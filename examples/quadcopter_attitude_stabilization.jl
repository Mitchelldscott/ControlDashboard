using Dash,
    DataFrames,
    PlotlyJS,
    StaticArrays,
    DifferentialEquations,
    ControlDashboard,
    ControlDashboard.ControlPanel

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
        integrator_error = SVector(0.0, 0.0, 0.0),
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
        new(
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
            integrator_error,
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

"""
    control_torque!(state, integral_error, kp, ki, kd) -> SVector{3,Float64}

Compute body-frame control torques `[τx, τy, τz]` for a quadcopter using a PID law
with proportional (`kp`), integral (`ki`), and derivative (`kd`) gains.

# Arguments

  - `state` :: `SVector{6,Float64}`
    The quadcopter attitude state vector:
    `["roll", "pitch", "yaw", "p", "q", "r"]`
    where `(roll, pitch, yaw)` are Euler angles [rad], and `(p, q, r)` are angular rates [rad/s].
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
    # state = ["roll", "pitch", "yaw", "p", "q", "r"]
    roll, pitch, yaw, p, q, r = state

    # Reference is zero, so error is just -state angles
    # Also pretends that there is no yaw sensor (usually there isn't)
    err = @SVector [-roll, -pitch, 0]

    # Derivative error = -angular rates (p,q,r)
    derr = @SVector [-p, -q, -r]

    # Update integral error
    integral_error = integral_error .+ err

    # PID control law
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

    # Build mixing matrix dynamically
    mixing = zeros(4, 4)
    for i in 1:4
        τ_i = cross(P[i], SVector(0, 0, kf))
        mixing[1, i] = kf
        mixing[2, i] = τ_i[1]
        mixing[3, i] = τ_i[2]
        mixing[4, i] = spin_dirs[i] * km
    end

    # Solve for squared speeds
    ω_sq = mixing \ SVector(thrust, τx, τy, τz)
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

    τ = SVector{3, Float64}(0.0, 0.0, 0.0)
    F_z = 0.0 # no gravity

    for i in 1:4
        ω² = w[i]^2
        Ti = kf * ω²
        Mi = km * ω² * spin_dirs[i]

        # Accumulate total thrust
        F_z += Ti
        τ += cross(motor_positions[i], SVector(0.0, 0.0, Ti)) + SVector(0.0, 0.0, Mi)
    end

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
    φ, θ, _ψ, p, q, r = u

    # Parameters extraction
    I_x, I_y, I_z = params.I_diag
    J_r = params.J_r
    A_r = params.Ar
    kp, ki, kd = params.kp, params.ki, params.kd
    integrator = params.integrator
    spin_dirs = params.spin_dirs

    # Delay controls (simulate startup)
    if t > 2.0
        τ_control = control_torque!(u, integrator, kp, ki, kd)
        ωs = motor_mixing(0.0, τ_control, params) # add throttle for altitude ctrl
        F, (τx, τy, τz) = aerodynamics(ωs, params)
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
    params = QuadcopterSimParameters(),
    template = "plotly_dark",
    frame_duration = 50,
)
    @assert all(name -> hasproperty(df, name), (:roll, :pitch, :yaw)) "DataFrame must contain :roll, :pitch, :yaw"
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
        QuadcopterSimParameters();
        shape = (3, 6),
        panel_style = Dict("display" => "flex", "alignItems" => "center"),
    )
end

"""
    initial_state(interfaces::Dict{String, Float64}) -> t_final, dt, x0, params, state_names

Produce the initial state of the quadcopter based on interface slider values.
Expected keys in `interfaces`: "roll", "pitch", "yaw",...
"""
function initialize_sim(inputs::NTuple{18, Any})
    names = (
        :t_final,
        :dt,
        :rpy_table,
        :pqr_table,
        :I_diag_table,
        :J_r,
        :Ar,
        :kp,
        :ki,
        :kd,
        :integrator_error_table,
        :m,
        :g,
        :L,
        :kf,
        :km,
        :motor_position_table,
        :spin_dir_table,
    )

    # Convert to named tuple for clarity
    kwargs = (; zip(names, inputs)...)

    # Fill in missing (nothing) entries
    # params_kw = merge(defaults, map(v -> isnothing(v[2]) ? (v[1] => get(defaults, v[1], v[2])) : v, pairs(kwargs)))

    params = QuadcopterSimParameters(; params_kw...)

    x0 = [
        deg2rad(params.rpy[1]),
        deg2rad(params.rpy[2]),
        deg2rad(params.rpy[3]),
        params.pqr[1],
        params.pqr[2],
        params.pqr[3],
    ]

    return params, x0
end

function quadcopter_simulation((params, x0))
    @info "Running Simulation"
    # run sim with RK4
    rk4_simulation(
        attitude_dynamics!,
        x0;
        t_final = params.t_final,
        dt = params.dt,
        params = params,
        state_names = ["roll", "pitch", "yaw", "p", "q", "r"],
    )
end

# --- Main execution ---
function main()
    run_dashboard(
        "Quadcopter Attitude Stabilizer",
        quadcopter_interfaces(),
        initialize_sim,
        quadcopter_simulation,
        Dict("main_view" => animate_quadcopter);
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
