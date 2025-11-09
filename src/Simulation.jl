module Simulation

using DataFrames, StaticArrays, OrdinaryDiffEq, SciMLBase

export sol_to_dataframe, rk4_simulation

"""
    sol_to_dataframe(sol::SciMLBase.AbstractODESolution; cols::AbstractVector)

Convert the solution of a `DifferentialEquations.jl` ODE problem into a `DataFrame`.

This function extracts the time vector and state trajectories from the ODE solution
and returns them in tabular form. It is primarily used to prepare simulation results
for visualization or further analysis (e.g., with Plotly or Dash).

# Arguments

  - `sol::SciMLBase.AbstractODESolution` : The solution returned by `DifferentialEquations.solve()`.
  - `cols::AbstractVector{<:AbstractString}` : Names for each state variable. The number of names must
    match the number of states in `sol`.

# Returns

  - `df::DataFrame` : A DataFrame containing the time vector and state trajectories.
    The first column is `:time`, followed by one column per state variable.

# Example

```julia
using DifferentialEquations, DataFrames

# Define and solve a simple ODE
f(u, p, t) = -0.5u
prob = ODEProblem(f, 1.0, (0.0, 5.0))
sol = solve(prob, Tsit5())

# Convert to DataFrame
df = sol_to_dataframe(sol; cols = ["u"])
first(df, 5)
```
"""
function sol_to_dataframe(sol::SciMLBase.AbstractODESolution; cols::AbstractVector)
    # Ensure solver finished successfully
    if !SciMLBase.successful_retcode(sol)
        error("Solution did not complete successfully. Return code: $(sol.retcode)")
    end

    # Convert solution to array form
    arr = Array(sol)  # each row = variable, each column = time sample
    nstates, nt = size(arr)

    # Check that names match state dimension
    if length(cols) != nstates
        error(
            "Number of columns ($cols) does not match the number of states ($nstates).",
        )
    end

    # Assemble DataFrame
    df = DataFrame(; time = sol.t)
    for (i, name) in enumerate(cols)
        df[!, Symbol(name)] = arr[i, :]
    end

    return df
end

"""
    rk4_simulation(system_dynamics, initial_state; t_final=5.0, dt=0.01, params=nothing, state_names=String[])

Simulate a system of ordinary differential equations (ODEs) over a fixed time horizon
using an explicit Runge-Kutta integrator (`Tsit5()`).

# Arguments

  - `system_dynamics`: A function `f!(du, u, p, t)` defining the system dynamics **in-place**.
  - `initial_state`: The initial state vector (or a function returning it).
  - `t_final`: Final simulation time (default = `5.0`).
  - `dt`: Time step used for saving results (default = `0.01`).
  - `params`: Optional parameters passed to the system dynamics (default = `nothing`).
  - `state_names`: Optional vector of state variable names for labeling output columns.

# Returns

  - `DataFrame` containing the simulation results with columns:

      + `time`: Simulation time points.
      + State variables (`x1, x2, ...` by default, or custom names if provided via `state_names`).

# Example

```julia
using DifferentialEquations, DataFrames

function pendulum!(du, u, p, t)
    θ, ω = u
    g, L = p
    du[1] = ω
    du[2] = -(g / L) * sin(θ)
end

u0 = [π / 2, 0.0]
p = (9.81, 1.0)

df = rk4_simulation(
    pendulum!,
    u0;
    t_final = 10.0,
    dt = 0.01,
    params = p,
    state_names = ["θ", "ω"],
)    # Define the time span and ODE problem
first(df, 5)
```
"""
function rk4_simulation(
    system_dynamics,
    initial_state;
    t_final::Real = 5.0,
    dt::Real = 0.01,
    params = nothing,
    state_names::AbstractVector{<:AbstractString} = String[],
)
    # Define the time span and ODE problem
    tspan = (zero(t_final), t_final)
    prob = SciMLBase.ODEProblem(system_dynamics, initial_state, tspan, params)

    # Integrate using a robust 5th-order explicit Runge–Kutta method
    sol = SciMLBase.solve(prob, OrdinaryDiffEq.Tsit5(); saveat = dt)

    # Convert to DataFrame
    return sol_to_dataframe(sol; cols = state_names)
end
end
