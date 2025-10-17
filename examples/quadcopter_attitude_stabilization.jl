using Dash,
    DataFrames,
    PlotlyJS,
    StaticArrays,
    DifferentialEquations

using ControlDashboard
using ControlDashboard.ControlPanel

struct QuadcopterSimParameters
    t_final::Float64             # Length of the simulation [s]
    dt::Float64                  # Sampling time of the simulation [s]
    rpy::SVector{3, Float64}      # Initial Attitude [degree]
    pqr::SVector{3, Float64}      # Initial Angular Rates [rad/s]
    I_diag::SVector{3, Float64}   # Diagonal inertia [kg·m^2]
    J_r::Float64                 # Rotor inertia [kg·m^2]
    Ar::Float64                  # Aerodynamic drag coefficient
    kp::Float64                  # Proportional gains
    ki::Float64                  # Integral gains
    kd::Float64                  # Derivative gains
    integrator::SVector{3, Float64} # Integral error state
    m::Float64                   # Mass [kg]
    g::Float64                   # Gravity [m/s^2]
    L::Float64                   # Arm length [m]
    kf::Float64                  # Thrust coefficient
    km::Float64                  # Drag torque coefficient
    motor_positions::SVector{4, SVector{3, Float64}}  # Motor positions in body frame
    spin_dirs::SVector{4, Int}           # spin directions: +1 for CCW, -1 for CW (used for yaw sign)

    function QuadcopterSimParameters(;
        t_final = 10.0,
        dt = 0.1,
        rpy = SVector(0.0, 0.0, 0.0),
        pqr = SVector(0.0, 0.0, 0.0),
        I_diag = SVector(1e-3, 1e-3, 1e-2),
        J_r = 6e-5,
        Ar = 1e-6,
        kp = 0.05,
        ki = 0.0,
        kd = 0.025,
        integrator = SVector(0.0, 0.0, 0.0),
        m = 0.05,
        g = 9.81,
        L = 0.1,
        kf = 1e-5,
        km = 2e-6,
        motor_positions = SVector(
            SVector(L / √2, L / √2, 0.0),   # M1: +x, +y
            SVector(L / √2, -L / √2, 0.0),   # M2: +x, -y
            SVector(-L / √2, -L / √2, 0.0),   # M3: -x, -y
            SVector(-L / √2, L / √2, 0.0),    # M4: -x, +y
        ),
        spin_dirs = SVector(1, -1, 1, -1),
    )
        return new(
            t_final,
            dt,
            rpy,
            pqr,
            I_diag,
            J_r,
            Ar,
            kp,
            ki,
            kd,
            integrator,
            m,
            g,
            L,
            kf,
            km,
            motor_positions,
            spin_dirs,
        )
    end
end

const SIM_SETTINGS = Ref{QuadcopterSimParameters}(QuadcopterSimParameters())

"""
    control_torque!(state, integral_error, kp, ki, kd) -> SVector{3,Float64}

Compute body-frame control torques `[τx, τy, τz]` for a quadcopter using a PID law
with proportional (`kp`), integral (`ki`), and derivative (`kd`) gains.

# Arguments

  - `state` :: `SVector{6,Float64}`
    The quadcopter attitude state vector:
    `["φ", "θ", "yaw", "p", "q", "r"]`
    where `(φ, θ, yaw)` are Euler angles [rad], and `(p, q, r)` are angular rates [rad/s].
  - `integral_error` :: `SVector{3,Float64}`
    Accumulated integral of attitude errors `[∫roll_err, ∫pitch_err, ∫yaw_err]`.
  - `kp`, `ki`, `kd` :: `SVector{3,Float64}` or scalar
    PID gains. Can be axis-specific vectors or scalars applied uniformly.

# Returns

  - `SVector{3,Float64}`
    Control torques `[τx, τy, τz]` in the body frame.

# Notes

  - The reference attitude is assumed to be **zero**, so the error is simply
    the negative of the current Euler angles, this function assumes only the
    roll and pitch measurements are available.
  - Derivative error is computed from the body angular rates `[-p, -q, -r]`.
  - The integral error is **updated in place** inside the function; if persistence
    across timesteps is required, the updated integral state must be stored externally.

# Role in Feedback Control Pipeline

(measured attitude) → [roll, pitch, yaw, p, q, r] → PID(state - reference) →
motor mixing → rotor thrusts → body torques [τx, τy, τz]

This function forms the **attitude controller** in the quadcopter feedback loop,    # state = ["roll", "pitch", "yaw", "p", "q", "r"]
stabilizing the vehicle around the desired orientation.
"""
function control_torque!(state, integral_error, kp, ki, kd)
    φ, θ, _, p, q, r = state
    err = @SVector [-φ, -θ, 0]
    derr = @SVector [-p, -q, -r]
    integral_error = integral_error .+ err
    return kp .* err + ki .* integral_error + kd .* derr
end

"""
    motor_mixing(thrust, τ_control, params) -> SVector{4, Float64}

Compute the squared motor speeds required to achieve a desired total thrust and body torques
for an arbitrary quadrotor geometry.

# Arguments

  - `thrust::Float64`: Desired total thrust (N) produced by all rotors combined.

  - `τ_control::SVector{3,Float64}`: Desired body torques `(τx, τy, τz)` in N·m.
  - `params`: Parameter struct containing:

      + `motor_positions::Vector{SVector{3,Float64}}`: Rotor position vectors in body frame (m).
      + `kf::Float64`: Thrust coefficient (N·s²/rad²).
      + `km::Float64`: Moment (drag) coefficient (N·m·s²/rad²).
      + `spin_dirs::NTuple{4,Float64}`: Spin direction of each rotor (+1 for CCW, −1 for CW).

# Returns

  - `SVector{4,Float64}`: The squared angular velocities `(ω₁², ω₂², ω₃², ω₄²)` of the motors
    required to generate the commanded thrust and torques. Negative values are clamped to zero.
"""
function motor_mixing(thrust, τ_control, params)
    τx, τy, τz = τ_control
    P = params.motor_positions
    kf, km = params.kf, params.km
    spin_dirs = params.spin_dirs

    τ = cross.(P, Ref(SVector(0.0, 0.0, kf)))   # τ_i = cross(P[i], [0,0,kf])
    mixing = @SMatrix [
        kf kf kf kf;
        τ[1][1] τ[2][1] τ[3][1] τ[4][1];
        τ[1][2] τ[2][2] τ[3][2] τ[4][2];
        spin_dirs[1]*km spin_dirs[2]*km spin_dirs[3]*km spin_dirs[4]*km
    ]

    # Solve for squared speeds
    ω_sq = pinv(mixing) * SVector{4}(thrust, τx, τy, τz)
    return map(x -> max(0.0, x), ω_sq)
end

"""
    aerodynamics(w_sq, params) -> SVector{3,Float64}

Compute the actual body-frame torques `[τx, τy, τz]` generated by the motors
spinning at given (squared) angular velocities.

This function models the physics of the quadcopter and is the "forward" model,
which is the inverse of the `motor_mixing` function. It's used within the
dynamics simulation to determine how the motor outputs affect the vehicle's state.

# Arguments

  - `w_sq` :: `NTuple{4, Float64}`
    A tuple of the squared angular velocities `(w1^2, w2^2, w3^2, w4^2)`.
  - `params` :: `QuadcopterSimParameters`
    A struct containing the quadcopter's physical parameters (`L`, `kf`, `km`).

# Returns

  - `SVector{3,Float64}`
    The resultant aerodynamic torques `[τx, τy, τz]` in the body frame.
"""
function aerodynamics(w, params)
    kf, km = params.kf, params.km
    motor_positions = params.motor_positions
    spin_dirs = params.spin_dirs

    F_z = 0.0 # no gravity
    w2 = w .^ 2  # element-wise square
    T = kf .* w2                     # thrust per motor
    M = km .* w2 .* spin_dirs        # moment per motor
    F_z = sum(T)                      # scalar
    τ = sum(cross.(motor_positions, @. SVector(0.0, 0.0, T)) .+ SVector.(0.0, 0.0, M))

    return F_z, τ
end

"""
    attitude_dynamics!(du, u, p, t)

Rigid-body rotational dynamics for a quadcopter.

  - u = [φ, θ, ψ, p, q, r] → roll, pitch, yaw and angular rates (Euler angles and Body Rates)
  - du = derivative
  - p = NamedTuple{(:torques, :I_diag, :J_r, :Ar, :Omega_r)} containing control moments and physical parameters.
"""
function attitude_dynamics!(du, u, params, t)
    # State extraction
    φ, θ, _, p, q, r = u

    # Parameters extraction
    I_x, I_y, I_z = params.I_diag
    J_r = params.J_r
    A_r = params.Ar
    kp, ki, kd = params.kp, params.ki, params.kd
    integrator = params.integrator
    spin_dirs = params.spin_dirs

    # Delay controls (simulate startup time)
    if t > 2.0
        τ_control = control_torque!(u, integrator, kp, ki, kd)
        ωs = motor_mixing(0.0, τ_control, params) # add throttle for altitude ctrl
        _, (τx, τy, τz) = aerodynamics(ωs, params)
    else
        ωs = (0.0, 0.0, 0.0, 0.0)
        (τx, τy, τz) = (0.0, 0.0, 0.0)
    end

    # --- 1. Attitude kinematics (Euler rates) ---
    sin_φ, cos_φ = sin(φ), cos(φ)
    cos_θ = cos(θ)

    # avoid singularity; if near gimbal-lock prefer quaternion model
    if abs(cos_θ) < 1e-6
        tan_θ = 0.0
        sec_θ = 0.0
    else
        tan_θ = tan(θ)
        sec_θ = 1.0 / cos_θ
    end

    du[1] = p + sin_φ * tan_θ * q + cos_φ * tan_θ * r
    du[2] = cos_φ * q - sin_φ * r
    du[3] = sin_φ * sec_θ * q + cos_φ * sec_θ * r

    # --- 2. Attitude dynamics (Newton-Euler) ---
    # Gyroscopic term: compute net rotor spin (signed)
    # sign convention: spin_dirs indicates +1 for CCW rotors (positive spin)
    Ω_net = sum(spin_dirs[i] * ωs[i] for i in 1:4)   # net rotor angular velocity (rad/s)

    # rate damping terms (Ar is linear damping of angular rates)
    du[4] = (I_y - I_z) / I_x * q * r - (J_r / I_x) * q * Ω_net + τx / I_x - A_r / I_x * p
    du[5] = (I_z - I_x) / I_y * p * r + (J_r / I_y) * p * Ω_net + τy / I_y - A_r / I_y * q
    du[6] = (I_x - I_y) / I_z * p * q + τz / I_z - A_r / I_z * r
end

"""
    relative_motor_positions(φ, θ, ψ, wing_length)

Compute the positions of the four quadcopter motors in 3D space, rotated according
to the quadcopter's roll (φ), pitch (θ), and yaw (ψ) angles. Positions are relative to the
center of mass (COM). Assumes '+' configuration.

# Arguments

  - `φ` : roll angle in radians
  - `θ` : pitch angle in radians
  - `ψ` : yaw angle in radians
  - `wing_length` : distance from COM to each motor (l)

# Returns

    # Rotation matrix R_IB (body → inertial) using ZYX (yaw-pitch-roll)

  - `Vector{Vector{Float64}}`: List of 4 motor positions `[x, y, z]` in the inertial frame after rotation.
"""
function relative_motor_positions(φ, θ, ψ, motor_positions)
    # Rotation matrix R_IB (body → inertial) using ZYX (yaw-pitch-roll)
    cψ, sψ = cos(ψ), sin(ψ)
    cθ, sθ = cos(θ), sin(θ)
    cφ, sφ = cos(φ), sin(φ)

    R_IB = @SMatrix [
        cθ*cψ sφ * sθ * cψ-cφ * sψ cφ * sθ * cψ+sφ * sψ;
        cθ*sψ sφ * sθ * sψ+cφ * cψ cφ * sθ * sψ-sφ * cψ;
        -sθ sφ*cθ cφ*cθ
    ]

    # Rotate body-frame motor positions into inertial frame
    P_inertial = [R_IB * P for P in motor_positions]
    # Return as a Vector of Float64 vectors
    return reduce(hcat, [collect(P) for P in P_inertial])
end

"""
    animate_quadcopter(df::DataFrame; params=QuadcopterSimParameters(), template="plotly_dark", frame_duration=50)

Build a PlotlyJS animation from a DataFrame `df` that must have columns:
:time, :roll, :pitch, :yaw

Each row produces a single frame; motor positions are computed by
`relative_motor_positions(roll, pitch, yaw, params)`. COM is at origin.

Returns a PlotlyJS.Plot ready to be used as the `figure` for `dcc_graph`.
"""
function animate_quadcopter(
    df;
    params = SIM_SETTINGS[],
    template = "plotly_dark",
    frame_duration = 50,
)
    wing_length = params.L
    quad_color = "black"
    prop_color = ["orange"; fill("red", 3)...]
    layout = Layout(;
        template = template,
        showlegend = false,
        scene = attr(;
            aspectmode = "manual",   # don't auto-stretch
            aspectratio = attr(; x = 1, y = 1, z = 1),  # equal scaling
            xaxis = attr(; range = [-3 * wing_length, 3 * wing_length]),
            yaxis = attr(; range = [-3 * wing_length, 3 * wing_length]),
            zaxis = attr(; range = [-3 * wing_length, 3 * wing_length]),
        ),
        updatemenus = [
            attr(;
                type = "buttons",
                showactive = false,
                buttons = [
                    attr(;
                        label = "Play",
                        method = "animate",
                        args = [
                            nothing,
                            attr(;
                                frame = attr(; duration = frame_duration, redraw = true),
                                fromcurrent = true,
                            ),
                        ],
                    ),
                ],
            ),
        ],
    )

    frames = PlotlyFrame[]

    # loop over timesteps in df
    for i in 1:nrow(df)
        roll = df.roll[i]
        pitch = df.pitch[i]
        yaw = df.yaw[i]

        # calculate motor positions from Euler angles
        motors = relative_motor_positions(roll, pitch, yaw, params.motor_positions)
        push!(
            frames,
            frame(;
                data = [
                    scatter3d(;
                        x = motors[1, :],
                        y = motors[2, :],
                        z = motors[3, :],
                        mode = "markers",
                        marker = attr(; size = 5),
                        line = attr(; width = 2, color = prop_color),
                        template = template,
                    ),
                    scatter3d(;
                        x = [motors[1, 1], motors[1, 3]],
                        y = [motors[2, 1], motors[2, 3]],
                        z = [motors[3, 1], motors[3, 3]],
                        mode = "lines",
                        line = attr(; width = 3, color = quad_color),
                        template = template,
                    ),
                    scatter3d(;
                        x = [motors[1, 2], motors[1, 4]],
                        y = [motors[2, 2], motors[2, 4]],
                        z = [motors[3, 2], motors[3, 4]],
                        mode = "lines",
                        line = attr(; width = 3, color = quad_color),
                        template = template,
                    ),
                ],
                layout = attr(;
                    annotations = [
                        attr(;
                            text = "Frame $i",
                            x = 0,
                            y = 1,  # top-left corner
                            xref = "paper",
                            yref = "paper",  # relative to entire figure
                            showarrow = false,
                            font = attr(; size = 16, color = "black"),
                        ),
                    ],
                ),
                name = "frame$i",
            ),
        )
    end

    # initial frame
    initial_trace = frames[1].data

    # Return the plot object (optional)
    return Plot(initial_trace, layout, frames)
end

# returns a Vector of components (same shape as your original quadcopter_interfaces)
function quadcopter_interfaces()
    return make_control_panel(
        SIM_SETTINGS[];
        component_style = Dict(
            "width" => "100%",
            "display" => "flex",
            "flex-direction" => "column",
            "align-items" => "stretch",
            "justify-content" => "center",
            "color" => "black",
            "padding" => "6px",
            "box-sizing" => "border-box",
        ),
        label_style = Dict(
            "font-weight" => "bold",
            "text-align" => "center",
            "margin-bottom" => "6px",
            "font-size" => "14px",
            "color" => "black",
            "display" => "block",
        ),
        panel_style = Dict(
            "display" => "grid",
            "grid-template-columns" => "repeat(6, 1fr)",  # 6 equal-width columns
            "grid-template-rows" => "repeat(3, auto)",    # 3 rows auto-sized
            "gap" => "16px",
            "width" => "100%",
            "height" => "auto",
            "padding" => "24px",
            "box-sizing" => "border-box",
            "background-color" => "#f8f8f8",
            "border-radius" => "12px",
            "box-shadow" => "0 4px 12px rgba(0, 0, 0, 0.1)",
            "justify-items" => "stretch",
            "align-items" => "start",
        ),
    )
end

"""
    quadcopter_simulation(interfaces::Dict{String, Float64}) -> t_final, dt, x0, params, state_names

Produce the initial state of the quadcopter based on interface slider values.
Expected keys in `interfaces`: "roll", "pitch", "yaw",...
"""
function quadcopter_simulation(inputs::NTuple{18, Any})
    names = fieldnames(QuadcopterSimParameters)

    # Convert to named tuple for clarity
    kwargs = (; zip(names, inputs)...)

    # For 1D vector table
    function read_num_vector(table)
        SVector{length(table)}([row[first(keys(row))] for row in table])
    end

    # For 2D table (matrix)
    function read_num_matrix(table)
        colnames = collect(keys(first(table)))
        nrows = length(table)
        ncols = length(colnames)
        SVector{nrows}(
            SVector{ncols}(row[name] for name in sort(colnames)) for
            row in table
        )
    end

    # Process each table into a simple SVector
    rpy = read_num_vector(kwargs.rpy)
    pqr = read_num_vector(kwargs.pqr)
    I_diag = read_num_vector(kwargs.I_diag)
    integrator = read_num_vector(kwargs.integrator)
    spin_dirs = read_num_vector(kwargs.spin_dirs)

    # Assuming motor_position contains vectors
    motor_positions = read_num_matrix(kwargs.motor_positions)
    # --- End Data Extraction ---

    # Construct the parameters struct using the corrected variable names and extracted data
    SIM_SETTINGS[] = QuadcopterSimParameters(;
        t_final = kwargs.t_final,
        dt = kwargs.dt,
        rpy = rpy,
        pqr = pqr,
        I_diag = I_diag,
        J_r = kwargs.J_r,
        Ar = kwargs.Ar,
        kp = kwargs.kp,
        ki = kwargs.ki,
        kd = kwargs.kd,
        integrator = integrator,
        m = kwargs.m,
        g = kwargs.g,
        L = kwargs.L,
        kf = kwargs.kf,
        km = kwargs.km,
        motor_positions = motor_positions,
        spin_dirs = spin_dirs,
    )

    # Initial state vector [phi, theta, psi, p, q, r]
    # Note: The constructor expects degrees for rpy, so we use the extracted values directly.
    # The deg2rad conversion is now done here for the initial state.
    x0 = [
        deg2rad(rpy[1]),
        deg2rad(rpy[2]),
        deg2rad(rpy[3]),
        pqr[1],
        pqr[2],
        pqr[3],
    ]

    @info "Running Simulation"

    # run sim with RK4
    rk4_simulation(
        attitude_dynamics!,
        x0;
        t_final = kwargs.t_final,
        dt = kwargs.dt,
        params = SIM_SETTINGS[],
        state_names = ["roll", "pitch", "yaw", "p", "q", "r"],
    )
end

# --- Main execution ---
function main()
    run_dashboard(
        "Quadcopter Attitude Stabilizer",
        quadcopter_interfaces(),
        quadcopter_simulation,
        Dict("main_view" => animate_quadcopter);
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
