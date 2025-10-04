using Dash, DataFrames, PlotlyJS, LinearAlgebra, DifferentialEquations, ControlDashboard, ControlDashboard.ControlPanel
using StaticArrays # Added for performance, standard practice in Julia rigid body dynamics [6, 7]

# Constants defined in the user's snippet
const I = Diagonal([0.01, 0.01, 0.02])  # kg·m²
const invI = inv(I)

# --- More Accurate Dynamics Constants (Derived from sources for typical quadrotor models) ---
# We define these constants globally for use in the simulation setup.
const J_ROTOR = 3.357e-5 # kg m^2 (Rotor inertia, Jr) [8]
const AR_DRAG = 0.0001 # Placeholder for aerodynamic resistance (Ar) [3]

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
    (τx, τy, τz) = control_torques!(u, params.integrator, params.kp, params.ki, params.kd)
    J_r = params.J_r
    A_r = params.Ar
    Omega_r = params.Omega_r

    # --- 1. Attitude Kinematics (Euler Angle Rates) ---
    # Transformation from body angular velocities (p, q, r) to Euler angle rates (φ̇, θ̇, ψ̇) [2]
    sin_φ, cos_φ = sin(φ), cos(φ)
    sin_θ, cos_θ = sin(θ), cos(θ)
    
    # Ensure cos_θ is non-zero (avoids gimbal lock/singularity near pitch θ = ±π/2)
    # If cos_θ is near zero, the model breaks down, requiring quaternion representation for accuracy [9].
    if abs(cos_θ) < 1e-6
        # Handle singularity condition (though this relies on the solver handling events/discontinuities 
        # which is a feature of DifferentialEquations.jl [10]) 
        # Setting rates to 0 or highly penalizing is common in numerical instability handling.
        # For typical flight simulation, we assume stability or use quaternions. We proceed assuming non-singular pitch.
        tan_θ = 0.0
        sec_θ = 0.0 
    else
        tan_θ = tan(θ)
        sec_θ = 1.0/cos_θ
    end

    du[1] = p + sin_φ * tan_θ * q + cos_φ * tan_θ * r  # φ̇ (Roll rate)
    du[2] = cos_φ * q - sin_φ * r                      # θ̇ (Pitch rate)
    du[3] = sin_φ * sec_θ * q + cos_φ * sec_θ * r      # ψ̇ (Yaw rate)

    # --- 2. Attitude Dynamics (Angular Accelerations, Newton-Euler Equations) ---
    # Based on Euler's equations for a rigid body, incorporating gyroscopic and aerodynamic drag terms [3]
    
    # Roll acceleration (ṗ)
    # ṗ = 1/Ix * [ (Iy - Iz)qr - Jr * q * Ωr + τx - Ar * p ] [120, Eq 2.11a, 55]
    du[4] = (I_y - I_z) / I_x * q * r - J_r / I_x * q * Omega_r + τx / I_x - A_r / I_x * p 

    # Pitch acceleration (q̇)
    # q̇ = 1/Iy * [ (Iz - Ix)pr + Jr * p * Ωr + τy - Ar * q ] [120, Eq 2.11b, 56]
    du[5] = (I_z - I_x) / I_y * p * r + J_r / I_y * p * Omega_r + τy / I_y - A_r / I_y * q

    # Yaw acceleration (ṙ)
    # ṙ = 1/Iz * [ (Ix - Iy)pq + τz - Ar * r ] [120, Eq 2.11c, 56]
    du[6] = (I_x - I_y) / I_z * p * q + τz / I_z - A_r / I_z * r
end

function control_torques!(state, integral_error, kp, ki, kd)
    # state = ["roll", "pitch", "yaw", "p", "q", "r"]
    roll, pitch, yaw, p, q, r = state

    # Reference is zero, so error is just -state angles
    err = @SVector [-roll, -pitch, -yaw]

    # Derivative error = -angular rates (p,q,r)
    derr = @SVector [-p, -q, -r]

    # Update integral error
    integral_error = integral_error .+ err

    # PID control law
    return kp .* err + ki .* integral_error + kd .* derr
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
- `Vector{Vector{Float64}}`: List of 4 motor positions `[x, y, z]` in the inertial frame after rotation.
"""
function relative_motor_positions(φ, θ, ψ, wing_length)
    l = wing_length

    # Positions in body frame ('+' configuration: +x, +y, -x, -y)
    P_body = [
        SVector(l, 0.0, 0.0),    # Front (x+)
        SVector(0.0, l, 0.0),    # Left (y+)
        SVector(-l, 0.0, 0.0),   # Rear (x-)
        SVector(0.0, -l, 0.0)    # Right (y-)
    ]

    # Rotation matrix R_IB (body → inertial) using ZYX (yaw-pitch-roll)
    cψ, sψ = cos(ψ), sin(ψ)
    cθ, sθ = cos(θ), sin(θ)
    cφ, sφ = cos(φ), sin(φ)

    R_IB = @SMatrix [
        cθ*cψ     sφ*sθ*cψ - cφ*sψ    cφ*sθ*cψ + sφ*sψ;
        cθ*sψ     sφ*sθ*sψ + cφ*cψ    cφ*sθ*sψ - sφ*cψ;
        -sθ       sφ*cθ               cφ*cθ
    ]

    # Rotate body-frame motor positions into inertial frame
    P_inertial = [R_IB * P for P in P_body]
    # Return as a Vector of Float64 vectors
    return reduce(hcat, [collect(P) for P in P_inertial])
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
    colorway = PlotlyJS.templates[template].layout[:colorway]
    quad_color = colorway[1]
    prop_color = colorway[2]
    layout = Layout(
        template=template,
        scene=attr(
            aspectmode="manual",   # don't auto-stretch
            aspectratio=attr(x=1, y=1, z=1),  # equal scaling
            xaxis=attr(range=[-3*wing_length, 3*wing_length]),
            yaxis=attr(range=[-3*wing_length, 3*wing_length]),
            zaxis=attr(range=[-3*wing_length, 3*wing_length])
        ),
        updatemenus=[attr(
            type="buttons",
            showactive=false,
            buttons=[attr(
                label="Play",
                method="animate",
                args=[nothing, attr(frame=attr(duration=frame_duration, redraw=true), fromcurrent=true)]
            )]
        )],
    )

    frames = PlotlyFrame[]

    # loop over timesteps in df
    for i in 1:nrow(df)
        roll  = df.roll[i]
        pitch = df.pitch[i]
        yaw   = df.yaw[i]

        # calculate motor positions from Euler angles
        motors = relative_motor_positions(roll, pitch, yaw, wing_length)
        push!(frames, frame(
            data=[
                scatter3d(
                    x=motors[1,:], 
                    y=motors[2,:], 
                    z=motors[3,:],
                    mode="markers",
                    marker=attr(size=5),
                    line=attr(width=2, color=prop_color),
                    template=template,
                ),
                scatter3d(
                    x=[motors[1,1], motors[1,3]],
                    y=[motors[2,1], motors[2,3]],
                    z=[motors[3,1], motors[3,3]],
                    mode="lines",
                    line=attr(width=3, color=quad_color),
                    template=template,
                ),
                scatter3d(
                    x=[motors[1,2], motors[1,4]],
                    y=[motors[2,2], motors[2,4]],
                    z=[motors[3,2], motors[3,4]],
                    mode="lines",
                    line=attr(width=3, color=quad_color),
                    template=template,
                )
            ],
            name="frame$i"
        ))
    end

    # initial frame
    initial_trace = frames[1].data

    return Plot(initial_trace, layout, frames)
end

# returns a Vector of components (same shape as your original quadcopter_interfaces)
function quadcopter_interfaces()
    return make_panel([
        Dict("component"=>"input", "label"=>"Duration", "id"=>"t_final", "value"=>10.0, "position"=>(1,1)),
        Dict("component"=>"input", "label"=>"Sample time", "id"=>"dt", "value"=>0.1, "position"=>(2,1)),
        Dict("component"=>"input", "label"=>"Roll", "id"=>"roll", "position"=>(1,2)),
        Dict("component"=>"input", "label"=>"Pitch", "id"=>"pitch", "position"=>(1,3)),
        Dict("component"=>"input", "label"=>"Yaw", "id"=>"yaw", "position"=>(1,4)),
        Dict("component"=>"input", "label"=>"P", "id"=>"p", "position"=>(2,2)),
        Dict("component"=>"input", "label"=>"Q", "id"=>"q", "position"=>(2,3)),
        Dict("component"=>"input", "label"=>"R", "id"=>"r", "position"=>(2,4)),
        Dict("component"=>"input", "label"=>"Kp", "id"=>"Kp", "value"=>0.1, "position"=>(3,1)),
        Dict("component"=>"input", "label"=>"Ki", "id"=>"Ki", "value"=>0.0, "position"=>(3,2)),
        Dict("component"=>"input", "label"=>"Kd", "id"=>"Kd", "value"=>0.05, "position"=>(3,3)),
    ]; shape=(3,4), panel_style=Dict(
        "display" => "flex",
        "alignItems" => "center"
    ))
end

"""
    initial_state(interfaces::Dict{String, Float64}) -> t_final, dt, x0, params, state_names

Produce the initial state of the quadcopter based on interface slider values.
Expected keys in `interfaces`: "roll", "pitch", "yaw",...
"""
function initial_state((t_final, dt, roll, pitch, yaw, p, q, r, kp, ki, kd))
    # Define the names of all the states to save into the df
    # Extract initial states from input
    state_names = ["roll", "pitch", "yaw", "p", "q", "r"]
    x0 = [roll, pitch, yaw, p, q, r]
    # Define parameters (p). Assuming zero constant control torques for an uncontrolled test 
    # and zero residual angular speed (Ωr) unless dynamically provided.
    params = (
        I_diag = diag(I),
        J_r = J_ROTOR,
        Ar = AR_DRAG,
        Omega_r = 0.0,
        kp = kp,
        ki = ki,
        kd = kd,
        integrator = SVector(0.0,0.0,0.0)
    )
    return t_final, dt, x0, params, state_names
end

function quadcopter_simulation((t_final, dt, x0, params, state_names))
    @info "Running Simulation"
    # run sim with RK4
    rk4_simulation(attitude_dynamics!, x0; t_final=t_final, dt=dt, params=params, state_names=state_names)
end

# --- Main execution ---
function main()
    app = initialize_dashboard("Quadcopter Attitude Controller"; interfaces=quadcopter_interfaces())
    set_callbacks!(
        app,
        initial_state, # Convert Control panel to initial state
        quadcopter_simulation, # Use RK4 to simulate
        Dict("main_view" => animate_quadcopter), # Render vizuals
        ["t_final", "dt", "roll", "pitch", "yaw", "p", "q", "r", "Kp", "Ki", "Kd"] # expected interfaces
    )
    run_server(app, "127.0.0.1", 8050)
end

main()

### References

```bibtex
@article{Rackauckas_Comparison_2018,
    author = {Rackauckas, Christopher},
    title = {A Comparison Between Differential Equation Solver Suites In MATLAB, R, Julia, Python, C, Mathematica, Maple, and Fortran},
    journal = {The Winnower},
    volume = {6},
    number = {e153459.98975},
    pages = {e153459.98975},
    year = {2018},
    doi = {10.15200/winn.153459.98975},
    note = {Cited in source}
}

@misc{JuliaDiscourse_RigidBodyDynamics_2020,
    author = {tkoolen and contributors},
    title = {Adding external force RigidBodyDynamics mechansms - General Usage},
    howpublished = {Julia Discourse},
    month = aug,
    year = {2020},
    note = {Cited in source}
}

@techreport{Determining_Products_Inertia,
    title = {Determining Products of Inertia for Small Scale UAVs},
    author = {{Undergraduate Student} and {Aerospace Engineer} and {Chief Scientist}},
    institution = {NASA Armstrong Flight Research Center},
    note = {Cited in source; specific authors and date not explicitly provided in excerpt.}
}

@misc{SciML_DifferentialEquations_jl_Docs,
    author = {{SciML}},
    title = {DifferentialEquations.jl: Multi-language suite for high-performance solvers of differential equations},
    howpublished = {Julia Packages Documentation},
    year = {2016--},
    note = {Ecosystem started in May 2016, cited in source}
}

@misc{JuliaDiscourse_DojoQuadcopter,
    title = {Dojo - simulating a quadcopter},
    howpublished = {Julia Discourse},
    note = {Cited in source}
}

@article{Idrissi_DynamicModelling_2020,
    author = {Idrissi, Moad and Annaz, Fawaz},
    title = {Dynamic Modelling and Analysis of a Quadrotor Based on Selected Physical Parameters},
    journal = {International Journal of Mechanical Engineering and Robotics Research},
    volume = {9},
    number = {6},
    month = jun,
    year = {2020},
    note = {Cited in source}
}

@misc{Seifner_FoundationInference_2024,
    author = {Seifner, Patrick and Cvejoski, Kostadin and Berghaus, David and Ojeda, C{\'e}sar and S{\'a}nchez, Rams{\'e}s J.},
    title = {Foundation Inference Models for Stochastic Differential Equations: A Transformer-based Approach for Zero-shot Function Estimation},
    archivePrefix = {arXiv},
    eprint = {2410.20587}, % Placeholder based on context, actual ID may differ
    year = {2024},
    note = {Under review; cited in source}
}

@online{DiffEq_GettingStarted,
    author = {{SciML} and {DifferentialEquations.jl contributors}},
    title = {Getting Started with Differential Equations in Julia},
    year = {2025},
    url = {https://docs.sciml.ai/DifferentialEquations/stable/tutorials/getting_started/}, % URL inferred from context
    note = {DifferentialEquations.jl documentation. Cited in source}
}

@online{SciML_GettingStarted_Overview,
    author = {{SciML}},
    title = {Overview of Julia's SciML: Quickly: What is Julia's SciML Ecosystem?},
    year = {2025},
    url = {https://docs.sciml.ai/SciML/stable/overview/}, % URL inferred from context
    note = {SciML documentation. Cited in source}
}

@online{JuliaHub_AerialVehicles,
    author = {{JuliaHub}},
    title = {Home $\cdot$ AerialVehicles.jl},
    url = {https://juliapackages.com/p/aerialvehicles},
    note = {Cited in source}
}

@online{JuliaHub_Multibody,
    author = {{JuliaHub} and {JuliaSim}},
    title = {Home $\cdot$ Multibody.jl},
    url = {https://juliapackages.com/p/multibody},
    note = {Cited in source}
}

@misc{rigidbodydynamicsjl,
    author = {Koolen, Twan and contributors},
    title = {{RigidBodyDynamics.jl}},
    year = {2016},
    url = {https://github.com/JuliaRobotics/RigidBodyDynamics.jl},
    note = {Explicitly cited format in source}
}

@misc{Suresh_HoveringControl,
    author = {Suresh, Harikrishnan and Sulficar, Abid and Desai, Vijay},
    title = {Hovering control of a quadcopter using linear and nonlinear techniques},
    note = {Cited in source}
}

@mastersthesis{Implementation_Quadcopter_Thesis,
    title = {Implementation and comparison of linearization-based and backstepping controllers for quadcopters},
    school = {Universitat Polit\`ecnica de Catalunya (UPCommons)},
    note = {Author name not explicitly provided in excerpt. Cited in source}
}

@inproceedings{Sells_JuliaBenchmark_FlightSim,
    author = {Sells, Ray},
    title = {Julia Programming Language Benchmark Using a Flight Simulation},
    booktitle = {Proceedings of the IEEE Aerospace Conference (Implied)},
    year = {TBD},
    organization = {NASA-MSFC},
    note = {Cited in source}
}

@techreport{Shen_Leok_LieGroupVI,
    author = {Shen, Xuefeng and Leok, Melvin},
    title = {LIE GROUP VARIATIONAL INTEGRATORS FOR RIGID BODY PROBLEMS USING QUATERNIONS},
    institution = {UCSD Math},
    note = {Cited in source}
}

@article{Liu_Modelling_PID_2024,
    author = {Liu, Yijin},
    title = {Modelling and simulation of quadcopter UAV based on PID control},
    journal = {Applied and Computational Engineering},
    volume = {75},
    pages = {202-210},
    month = jul,
    year = {2024},
    doi = {10.54254/2755-2721/75/20240539},
    note = {Cited in source}
}

@phdthesis{Nishimura_OnlineTrajectory_2021,
    author = {Nishimura, Haruki},
    title = {ONLINE TRAJECTORY PLANNING ALGORITHMS FOR ROBOTIC SYSTEMS UNDER UNCERTAINTY IN INTERACTIVE ENVIRONMENTS},
    school = {Stanford University},
    month = aug,
    year = {2021},
    note = {Cited in source}
}

@techreport{Wolter_Quadcopter_Report_2015,
    author = {Wolter, Moritz (Group 19)},
    title = {Quadcopter exercise - Simulation Report},
    institution = {Unspecified},
    month = jan,
    year = {2015},
    note = {Cited in source}
}

@techreport{Quadrotor_Modeling_Control_Unspecified,
    title = {Quadrotor Modeling and Control},
    institution = {Space Imaging and Optical Systems Lab},
    note = {Cited in source}
}

@online{RBD_QuickStartGuide,
    author = {{JuliaRobotics}},
    title = {Quick start guide $\cdot$ RigidBodyDynamics.jl},
    url = {https://juliarobotics.github.io/RigidBodyDynamics.jl/stable/quickstart/}, % URL inferred from context
    note = {RigidBodyDynamics.jl Documentation. Cited in source}
}

@misc{Review_PID_UAVs,
    title = {Review of PID Controller Applications for UAVs},
    note = {Cited in source}
}

@misc{JuliaDiscourse_StiffSolvers_2022,
    author = {graeffd and baggepinnen and isaacsas},
    title = {Solver for stiff differential equations},
    howpublished = {Julia Programming Language Discourse},
    month = may,
    year = {2022},
    note = {Cited in source}
}

@online{Ghosh_Solving_DE_Julia_2025,
    author = {Ghosh, Ritobrata},
    title = {Solving First Order Differential Equations with Julia},
    month = {mar},
    year = {2025},
    url = {https://ritog.github.io/posts/1st-order-DE-julia/1st_order_DE_julia.html},
    note = {Cited in source}
}

@online{AerialVehicles_Stabilizing,
    author = {{AerialVehicles.jl contributors}},
    title = {Stabilizing a falling quadrotor $\cdot$ AerialVehicles.jl},
    month = sep,
    year = {2023},
    note = {Cited in source}
}

@online{RoboticsSE_EulerAngles,
    title = {Treatment of euler angles in quadcopter control},
    howpublished = {Robotics Stack Exchange},
    note = {Cited in source}
}

@misc{Snyder_JuliaCon2020_DroneAutonomy,
    author = {Snyder, Kerry},
    title = {Rapid Commercialization of Drone Autonomy using Julia},
    booktitle = {JuliaCon 2020},
    howpublished = {Video presentation},
    year = {2020},
    note = {Cited in source}
}

@misc{Ferrolho_JuliaCon2023_LeggedRobots,
    author = {Ferrolho, Henrique},
    title = {Using Julia to Optimise Trajectories for Robots with Legs},
    booktitle = {JuliaCon 2023},
    howpublished = {Video presentation},
    year = {2023},
    note = {Cited in source}
}

@inproceedings{Ferrolho_InverseDynamics_2021,
    author = {Ferrolho, Henrique and Ivan, Vladimir and Merkt, Wolfgang and Havoutis, Ioannis and Vijayakumar, Sethu},
    title = {Inverse Dynamics vs. Forward Dynamics in Direct Transcription Formulations for Trajectory Optimization},
    booktitle = {2021 IEEE International Conference on Robotics and Automation (ICRA)},
    year = {2021},
    eprint = {2010.05359},
    archivePrefix = {arXiv},
    note = {Cited in source}
}

@misc{Saroufim_RigidBodyDynamicsTutorial,
    author = {Saroufim, msaroufim},
    title = {RigidBodyDynamicsTutorial: A tutorial on Rigid Body Dynamics in Julia},
    howpublished = {GitHub Repository/Notes},
    note = {Cited in source}
}

@misc{Plotly_Dash_jl_GitHub,
    author = {{plotly}},
    title = {Dash.jl: Dash for Julia - A Julia interface to the Dash ecosystem},
    howpublished = {GitHub Repository},
    note = {Cited in source}
}
```