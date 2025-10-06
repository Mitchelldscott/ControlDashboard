# ControlDashboard.jl
# A generalized Julia module to create an interactive control simulation dashboard using Dash.jl.
module ControlDashboard

# Import required libraries
using Dash, DataFrames, StaticArrays, DifferentialEquations

include("ControlPanel.jl")

# Export the main dashboard function so it can be used by other modules.
export initialize_dashboard, set_callbacks!, sol_to_dataframe, rk4_simulation, run_dashboard

"""
    initialize_dashboard(title; interfaces=make_interfaces(), graphs=make_graphs()) -> DashApp

Factory function that initializes a Dash app with the provided styles, interfaces and graphs.
"""
function initialize_dashboard(
    title;
    title_style = Dict("textAlign" => "center"),
    control_panel = sample_time_and_duration_sliders(),
    graphs = [dcc_graph(id = "main_view")],
    external_stylesheets = ["https://bootswatch.com/5/darkly/bootstrap.min.css"],
)
    app = dash(external_stylesheets = external_stylesheets)

    app.layout = html_div() do
        html_h1(title, style = title_style),
        control_panel,  # Control Panel
        graphs...    # vizualizations
    end

    return app
end

"""
    control_dashboard(state_factory, run_simulation, make_figure)

Register callbacks for interactive components of a generic Dash app.

This function uses a set of input functions to initialize and run a 
simulation, 

# Arguments
- `state_factory::Function`: A function that takes no arguments and returns a fresh 
initial state object (e.g., a Dict or a custom struct with initial conditions).
- `run_simulation::Function`: A function that executes the simulation. It should accept
the initial state as input and return a dataframe.
- `make_figure::Function`: A function that takes the simulation data and returns a
figure to display.
"""
function set_callbacks!(app, state_factory, run_simulation, figures, interfaces)
    # Define the callback function that links the sliders to the graph.
    # When a slider value changes, this function is triggered.
    for (figure_name, renderer) in figures
        callback!(
            app,
            Output(figure_name, "figure"),
            [Input(interface, "value") for interface in interfaces],
        ) do interfaces...
            df = run_simulation(state_factory(interfaces))
            return renderer(df)
        end
    end
end

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
function sol_to_dataframe(sol; cols = nothing)
    arr = Array(sol)  # each column is a state, rows = time steps
    df = DataFrame(time = sol.t) # sol should always have t
    # Auto-generate names if not provided
    if isnothing(cols)
        cols = ["x$(i)" for i = 1:size(arr, 1)]
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
    - `time` = simulation time points  
    - state columns extracted from the solution (`x1, x2, ...` by default, or renamed if using `sol_to_dataframe` with labels).
"""
function rk4_simulation(
    system_dynamics,
    initial_state;
    t_final = 5.0,
    dt = 0.01,
    params = nothing,
    state_names = nothing,
)
    # Use SVector for initial state for performance optimization recommended in Julia ODE/Robotics ecosystems [6]
    tspan = (0.0, t_final)
    # Define the ODE problem
    prob = ODEProblem(system_dynamics, initial_state, tspan)
    # Solve the ODE. Using Tsit5(), a powerful explicit Runge-Kutta method often effective 
    # for non-stiff dynamics like this attitude model [12, 16].
    sol = solve(prob, Tsit5(), saveat = dt, p = params)
    # Process results into DataFrame
    return sol_to_dataframe(sol; cols = state_names)
end

"""
    run_dashboard(control_panel::Vector, initialize_sim::Function, simulate::Function,
                  renderers::Dict{String, Function};
                  host::String="127.0.0.1", port::Int=8050, title::String="Control Dashboard")

Initialize and run an interactive dashboard from a pre-built control panel.

This function provides a high-level entry point for launching a `ControlDashboard`
application. The user supplies a fully constructed control panel (typically
from [`make_control_panel`](@ref)), simulation callbacks, and rendering
functions. The dashboard is assembled, callbacks are registered, and a local
server is started.

# Arguments
- `control_panel`: A `Vector` of Dash HTML or core components representing the UI
  layout (e.g. created by [`make_control_panel`](@ref)).
- `initialize_sim`: A function that converts control panel inputs into initial
  simulation states and parameter dictionaries.
- `simulate`: The main simulation function advancing system dynamics.
- `renderers`: A `Dict{String, Function}` mapping component IDs to functions that
  render plots, animations, or visual outputs.
- `host`: The server host address (default `"127.0.0.1"`).
- `port`: The server port (default `8050`).
- `title`: The title shown in the browser tab and top of the dashboard (default `"Control Dashboard"`).

# Returns
- A `Dash.DashApp` instance representing the running dashboard.

# Example
```julia
panel = make_control_panel(QuadCopterParams())
app = run_dashboard(
    panel,
    initialize_sim,
    quadcopter_simulation,
    Dict("main_view" => animate_quadcopter);
    title = "Quadcopter Attitude Stabilizer",
)
'''
This sets up a "Quadcopter Attitude Stabilizer" dashboard, runs the control
loop, and serves it on http://127.0.0.1:8050.
"""
function run_dashboard(
    control_panel::Vector,
    initialize_sim::Function,
    simulate::Function,
    renderers::Dict{String, Function};
    host::String = "127.0.0.1",
    port::Int = 8050,
    title::String = "Control Dashboard"
)
    # Step 1: Initialize dashboard app with user-specified interfaces
    app = initialize_dashboard(title; control_panel = control_panel, graphs = [dcc_graph(id = "main_view")])
    
    # Step 2: Register simulation and renderer callbacks
    set_callbacks!(
        app,
        initialize_sim,
        simulate,
        renderers,
        get_component_ids(control_panel),
    )

    # Step 4: Run the dashboard server
    run_server(app, host, port)

    return nothing
end
end # module ControlDashboards
