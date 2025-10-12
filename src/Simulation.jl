module Simulation

using DataFrames, StaticArrays, DifferentialEquations

export sol_to_dataframe, rk4_simulation

"""
    sol_to_dataframe(sol; state_names)

Construct a DataFrame from the solution of a DifferentialEquations.ODEProblem.

This is a helper function for simulations that use DifferentialEquations.jl.

# Arguments

  - `sol::` : The output of `DifferentialEquations.solve()`.
  - `cols::Vector{String}` : Column names of the DataFrame.

# Returns

  - `df::DataFrame` : States and sample times extracted from a sol.
"""
function sol_to_dataframe(sol; cols::AbstractVector = string[])
    arr = Array(sol)  # each column is a state, rows = time steps
    df = DataFrame(; time = sol.t) # sol should always have t
    # Auto-generate names if not provided
    if isnothing(cols)
        cols = ["x$(i)" for i in 1:size(arr, 1)]
    end
    for (i, name) in enumerate(cols)
        df[!, Symbol(name)] = arr[i, :]  # grab the i-th state trajectory
    end
    return df
end

"""
    rk4_simulation(system_dynamics, initial_state; t_final=5.0, dt=0.01, p=nothing)

Simulate a system of ordinary differential equations (ODEs) over a given time horizon.

# Arguments

  - `system_dynamics`: A function `f!(du, u, p, t)` defining the system dynamics in-place.
  - `initial_state`: Vector-like object or function that prepares simulation data from a tuple of interfaces.
  - `t_final`: (default = 10.0) Final simulation time.
  - `dt`: (default = 0.01) Time step for saving results.
  - `p`: (default = nothing) Optional parameters to pass to the dynamics.

# Returns

  - `DataFrame` containing simulation results with columns:

      + `time` = simulation time points
      + state columns extracted from the solution (`x1, x2, ...` by default, or renamed if using `sol_to_dataframe` with labels).
"""
function rk4_simulation(
    system_dynamics,
    initial_state;
    t_final::Number = 5.0,
    dt::Number = 0.01,
    params = nothing,
    state_names::AbstractVector = String[],
)
    # Use SVector for initial state for performance optimization recommended in Julia ODE/Robotics ecosystems [6]
    tspan = (0.0, t_final)
    # Define the ODE problem
    prob = ODEProblem(system_dynamics, initial_state, tspan)
    # Solve the ODE. Using Tsit5(), a powerful explicit Runge-Kutta method often effective 
    # for non-stiff dynamics like this attitude model [12, 16].
    sol = solve(prob, Tsit5(); saveat = dt, p = params)
    # Process results into DataFrame
    return sol_to_dataframe(sol; cols = state_names)
end
end # module Simulation
