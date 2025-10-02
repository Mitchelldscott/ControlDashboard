using Dash, DataFrames, PlotlyJS, LinearAlgebra, DifferentialEquations, ControlDashboard, ControlDashboard.ControlPanel

# Quadcopter inertia (example values)
const I = Diagonal([0.01, 0.01, 0.02])  # kg·m²
const invI = inv(I)

"""
    initial_state(interfaces::Dict{String, Float64}) -> Dict{String, Float64}

Produce the initial state of the quadcopter based on interface slider values.
Expected keys in `interfaces`: "roll", "pitch", "yaw", "motor_speeds" (array of 4 floats)
"""
function initial_state((roll, pitch, yaw))
    return Dict(
        "roll" => roll,
        "pitch" => pitch,
        "yaw" => yaw
    )
end

"""
    attitude_dynamics!(du, u, p, t)

Rigid-body rotational dynamics for a quadcopter.
- u = [φ, θ, ψ, p, q, r] → roll, pitch, yaw and angular rates
- du = derivative
- p = Dict("torques" => [τx, τy, τz])
"""
function attitude_dynamics!(du, u, p, _t)
    φ, θ, _ψ, p_rate, q_rate, r_rate = u
    τ = zeros(3)

    # Euler angles rates
    du[1] = p_rate + q_rate * sin(φ) * tan(θ) + r_rate * cos(φ) * tan(θ)
    du[2] = q_rate * cos(φ) - r_rate * sin(φ)
    du[3] = q_rate * sin(φ)/cos(θ) + r_rate * cos(φ)/cos(θ)

    # Angular acceleration (rigid body)
    omega = [p_rate, q_rate, r_rate]
    domega = invI * (τ - cross(omega, I*omega))

    du[4:6] .= domega
end

"""
    quadcopter_simulation(state; t_final=5.0, dt=0.01)

Simulate the attitude of a quadcopter given initial roll, pitch, yaw and angular rates.
Returns a DataFrame with columns: :time, :roll, :pitch, :yaw
"""
function quadcopter_simulation(state; t_final=5.0, dt=0.01)
    # Initial state: [roll, pitch, yaw, p, q, r]
    u0 = [
        state["roll"],
        state["pitch"],
        state["yaw"],
        0.0, 0.0, 0.0  # initial angular rates
    ]

    tspan = (0.0, t_final)
    prob = ODEProblem(attitude_dynamics!, u0, tspan)
    sol = solve(prob, Tsit5(), dt=dt)  # Tsit5 = 5th order Runge-Kutta

    return DataFrame(
        time = sol.t,
        roll = sol[1,:],
        pitch = sol[2,:],
        yaw = sol[3,:]
    )
end

"""
    relative_motor_positions(initial_conditions::Dict{String, Float64}) -> Vector{Vector{Float64}}

Compute the positions of the four quadcopter motors in 3D space, rotated according
to the quadcopter's roll, pitch, and yaw angles. Positions are relative to the
center of mass (COM) and assumed to be xy-symetric.

# Arguments
- `initial_conditions::Dict{String, Any}`: Dictionary with keys:
    - `φ` : roll angle in radians
    - `θ` : pitch angle in radians
    - `ψ` : yaw angle in radians
    - `wing_length` : distance from COM to each motor
# Returns
- `Vector{Vector{Float64}}`: List of 4 motor positions `[x, y, z]` after rotation.
"""
function relative_motor_positions(φ, θ, ψ, wing_length)
    # Rotation matrices (ZYX convention: yaw-pitch-roll)
    R_roll = [1 0 0; 0 cos(φ) -sin(φ); 0 sin(φ) cos(φ)]
    R_pitch = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
    R_yaw = [cos(ψ) -sin(ψ) 0; sin(ψ) cos(ψ) 0; 0 0 1]

    # Combined rotation
    R = R_yaw * R_pitch * R_roll

    # Motor positions in body frame (X forward)
    motor_positions = [
        [ 1.0,  1.0, 0.0] .* wing_length, # Front left 
        [ 1.0, -1.0, 0.0] .* wing_length, # Front right
        [-1.0, -1.0, 0.0] .* wing_length, # Rear right
        [-1.0,  1.0, 0.0] .* wing_length, # Rear left
    ]

    # Apply rotation
    rotated_positions = [R * pos for pos in motor_positions]

    # return as 4×3 matrix
    return reduce(vcat, (p' for p in rotated_positions))
end

"""
    animate_quadcopter(df::DataFrame; wing_length=1.0, template="plotly_dark", frame_duration=50)

Build a PlotlyJS animation from a DataFrame `df` that must have columns:
  :time, :roll, :pitch, :yaw

Each row produces a single frame; motor positions are computed by
`relative_motor_positions(roll, pitch, yaw, wing_length)`. COM is at origin.

Returns a PlotlyJS.Plot ready to be used as the `figure` for `dcc_graph`.
"""
function animate_quadcopter(df; wing_length=0.1, template="plotly_dark", frame_duration=50)

    @assert all(name -> hasproperty(df, name), (:roll, :pitch, :yaw)) "DataFrame must contain :roll, :pitch, :yaw"

    frames = PlotlyFrame[]

    # loop over timesteps in df
    for i in 1:nrow(df)
        roll  = df.roll[i]
        pitch = df.pitch[i]
        yaw   = df.yaw[i]

        # calculate motor positions from Euler angles
        motors = relative_motor_positions(roll, pitch, yaw, wing_length)
        
        push!(frames, frame(
            data=[scatter3d(
                x=motors[:,1], 
                y=motors[:,2], 
                z=motors[:,3],
                mode="markers",
                marker=attr(size=5, color="blue"),
                line=attr(width=2)
            )],
            name="frame$i"
        ))
    end

    # initial frame
    initial_trace = frames[1].data

    layout = Layout(
        updatemenus=[attr(
            type="buttons",
            showactive=false,
            buttons=[attr(
                label="Play",
                method="animate",
                args=[nothing, attr(frame=attr(duration=frame_duration, redraw=true), fromcurrent=true)]
            )]
        )],
        template=template,
    )

    return Plot(initial_trace, layout, frames)
end

# returns a Vector of components (same shape as your original quadcopter_interfaces)
function quadcopter_interfaces()
    return make_panel([
        Dict("component"=>"input", "label"=>"Roll", "id"=>"roll", "position"=>(1,1)),
        Dict("component"=>"input", "label"=>"Pitch", "id"=>"pitch", "position"=>(1,2)),
        Dict("component"=>"input", "label"=>"Yaw", "id"=>"yaw", "position"=>(1,3)),
        Dict("component"=>"input", "label"=>"P", "id"=>"p", "position"=>(2,1)),
        Dict("component"=>"input", "label"=>"Q", "id"=>"q", "position"=>(2,2)),
        Dict("component"=>"input", "label"=>"R", "id"=>"r", "position"=>(2,3)),
    ]; shape=(2,3), row_style=Dict(
        "display" => "flex",
        "alignItems" => "center",
        "gap" => "10px",
        "marginBottom" => "8px"
    ))
end

# --- Main execution ---
function main()
    app = initialize_dashboard("Quadcopter Attitude Controller"; interfaces=quadcopter_interfaces())
    app = set_callbacks(
        app,
        initial_state,
        quadcopter_simulation,
        Dict("main_view" => animate_quadcopter),
        ["roll", "pitch", "yaw"]
    )
    run_server(app, "127.0.0.1", 8050)
end

main()
