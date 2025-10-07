# ControlDashboard.jl
# A generalized Julia module to create an interactive control simulation dashboard using Dash.jl.
module ControlDashboard

# Import required libraries
using Dash, DataFrames, StaticArrays, DifferentialEquations

# Sub-modules
include("ControlPanel.jl")
using .ControlPanel: make_panel, make_control_panel, get_component_ids

include("Simulation.jl")
using .Simulation: sol_to_dataframe, rk4_simulation


# Export the main dashboard function so it can be used by other modules.
export 
    initialize_dashboard, set_callbacks!, run_dashboard,
    make_panel, make_control_panel, get_component_ids,
    sol_to_dataframe, rk4_simulation

"""
    initialize_dashboard(title; interfaces=make_interfaces(), graphs=make_graphs()) -> DashApp

Factory function that initializes a Dash app with the provided styles, interfaces and graphs.
"""
function initialize_dashboard(
    title::AbstractString;
    title_style::Dict = Dict("textAlign" => "center"),
    control_panel::AbstractVector = sample_time_and_duration_sliders(),
    views::AbstractVector = [dcc_graph(id = "main_view")],
    external_stylesheets::AbstractVector = ["https://bootswatch.com/5/darkly/bootstrap.min.css"],
)
    app = dash(external_stylesheets = external_stylesheets)

    app.layout = html_div() do
        html_h1(title, style = title_style),
        control_panel,  # Control Panel
        views...    # vizualizations
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
function set_callbacks!(
    app, 
    state_factory::Function, 
    run_simulation::Function, 
    figures::AbstractDict, 
    interfaces::AbstractVector
)
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
    run_dashboard(title::string, control_panel::Vector, initialize_sim::Function, 
                  simulate::Function, renderers::Dict{String, Function};
                  host::String="127.0.0.1", port::Int=8050, views=[dcc_graph(id="main_view")])

Initialize and run an interactive dashboard from a pre-built control panel.

This function provides a high-level entry point for launching a `ControlDashboard`
application. The user supplies a fully constructed control panel (typically
from [`make_control_panel`](@ref)), simulation callbacks, and rendering
functions. The dashboard is assembled, callbacks are registered, and a local
server is started.

# Arguments
- `title`: The title shown in the browser tab and top of the dashboard.
- `control_panel`: A `Vector` of Dash HTML or core components representing the UI
  layout (e.g. created by [`make_control_panel`](@ref)).
- `initialize_sim`: A function that converts control panel inputs into initial
  simulation states and parameter dictionaries.
- `simulate`: The main simulation function advancing system dynamics.
- `renderers`: A `Dict{String, Function}` mapping component IDs to functions that
  render plots, animations, or visual outputs.
- `host`: The server host address (default `"127.0.0.1"`).
- `port`: The server port (default `8050`).
- `views`: The figures that callbacks will render to.
# Returns
- A `Dash.DashApp` instance representing the running dashboard.

# Example
```julia
panel = make_control_panel(QuadCopterParams())
app = run_dashboard(
    "Quadcopter Attitude Stabilizer",
    panel,
    initialize_sim,
    quadcopter_simulation,
    Dict("main_view" => animate_quadcopter);
)
'''
This sets up a "Quadcopter Attitude Stabilizer" dashboard, runs the control
loop, and serves it on http://127.0.0.1:8050.
"""
function run_dashboard(
    title::AbstractString,
    control_panel::AbstractVector,
    initialize_sim::Function,
    simulate::Function,
    renderers::AbstractDict;
    host::AbstractString = "127.0.0.1",
    port::Integer = 8050,
    views::AbstractVector = [dcc_graph(id = "main_view")]
)
    # Step 1: Initialize dashboard app with user-specified interfaces
    app = initialize_dashboard(title; control_panel = control_panel, views = views)
    
    # Step 2: Register simulation and renderer callbacks
    set_callbacks!(
        app,
        initialize_sim,
        simulate,
        renderers,
        get_component_ids(control_panel)
    )

    # Step 3: Run the dashboard server
    run_server(app, host, port)

    return nothing
end
end # module ControlDashboard
