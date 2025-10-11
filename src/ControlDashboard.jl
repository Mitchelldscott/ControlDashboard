# ControlDashboard.jl
# A generalized Julia module to create an interactive control simulation dashboard using Dash.jl.
module ControlDashboard

# Import required libraries
using Dash, DataFrames, StaticArrays, DifferentialEquations

# Sub-modules
include("ControlPanel.jl")
using .ControlPanel: make_panel, make_control_panel, get_interactive_components

include("Simulation.jl")
using .Simulation: sol_to_dataframe, rk4_simulation

export 
    initialize_dashboard, set_callbacks!, run_dashboard,
    make_panel, make_control_panel, get_interactive_components,
    sol_to_dataframe, rk4_simulation

"""
    initialize_dashboard(title; title_style, control_panel, views, external_stylesheets)

Create and initialize a Dash application layout with a title, control panel, and visualization views.

# Arguments
- `title`: Page title displayed at the top of the dashboard.
- `title_style`: CSS style dictionary for the title (default centers text).
- `control_panel`: Vector of Dash components defining the interactive control panel.
- `views`: Vector of Dash visualization components (e.g., graphs, tables).
- `external_stylesheets`: URLs of external CSS stylesheets applied to the dashboard.

# Returns
- A configured `Dash` app instance with the specified layout.
"""
function initialize_dashboard(title::AbstractString;
    title_style::Dict = Dict("textAlign" => "center"),
    control_panel::AbstractVector = sample_time_and_duration_sliders(),
    views::AbstractVector = [dcc_graph(id = "main_view")],
    external_stylesheets::AbstractVector = ["https://bootswatch.com/5/darkly/bootstrap.min.css"])

    app = dash(external_stylesheets = external_stylesheets)
    
    app.layout = html_div() do
        html_h1(title, style = title_style),
        control_panel,  # Control Panel
        views...    # vizualizations
    end

    return app
end

"""
    set_callbacks!(app, state_factory, run_simulation, figures, interfaces)

Register reactive Dash callbacks linking UI components to visualization updates.

# Arguments
- `app`: The Dash application instance where callbacks will be registered.
- `state_factory`: Function that constructs a simulation or model state from user interface inputs.
- `run_simulation`: Function that runs the simulation or computation given the constructed state.
- `figures`: Dictionary mapping figure IDs to renderer functions that convert simulation results (e.g., `DataFrame`) into plotly figures.
- `interfaces`: Vector of `(id, field)` pairs identifying interactive components and the fields to observe (e.g., `("slider_1", "value")`).

# Behavior
For each figure defined in `figures`, a Dash callback is created.  
Whenever any input in `interfaces` changes, the callback:
1. Builds the state via `state_factory(interfaces)`.
2. Executes `run_simulation` on that state.
3. Updates the corresponding figure using its renderer.

If `figures` or `interfaces` are empty, no callbacks are registered.
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
    if length(interfaces) === 0 return end
    for (figure_name, renderer) in figures
        callback!(
            app,
            Output(figure_name, "figure"),
            [Input(name, field) for (name, field) in interfaces],
        ) do interfaces...
            df = run_simulation(state_factory(interfaces))
            return renderer(df)
        end
    end
end

"""
    run_dashboard(title, control_panel, initialize_sim, simulate, renderers;
                  host = "127.0.0.1", port = 8050,
                  views = [dcc_graph(id = "main_view")],
                  external_stylesheets = ["https://bootswatch.com/5/darkly/bootstrap.min.css"])

Launch an interactive Dash-based simulation dashboard that connects UI controls to visualization updates.

# Arguments
- `title`: Title displayed at the top of the dashboard.
- `control_panel`: Vector of Dash components defining the interactive control panel.
- `initialize_sim`: Function that constructs an initial simulation state from UI inputs.
- `simulate`: Function that runs the simulation or computation given a constructed state.
- `renderers`: Dictionary mapping figure IDs to rendering functions that convert simulation outputs (e.g., data tables or DataFrames) into Plotly figures.

# Keyword Arguments
- `host`: Host address for the dashboard server (default `"127.0.0.1"`).
- `port`: Port number for the dashboard server (default `8050`).
- `views`: Vector of Dash visualization components to display (default a single graph view).
- `external_stylesheets`: URLs of external CSS stylesheets applied to the dashboard theme.

# Behavior
1. Initializes the Dash app layout using `initialize_dashboard`.
2. Registers interactive callbacks that connect control panel inputs to simulation and visualization updates.
3. Starts the local web server and serves the dashboard.

# Returns
Nothing. Runs the dashboard server in-place.
"""
function run_dashboard(title::AbstractString,
    control_panel::AbstractVector,
    initialize_sim::Function,
    simulate::Function,
    renderers::AbstractDict;
    host::AbstractString = "127.0.0.1",
    port::Integer = 8050,
    views::AbstractVector = [dcc_graph(id = "main_view")],
    external_stylesheets::AbstractVector = ["https://bootswatch.com/5/darkly/bootstrap.min.css"])

    # Step 1: Initialize dashboard app with user-specified interfaces
    app = initialize_dashboard(title; 
        control_panel = control_panel, 
        views = views, 
        external_stylesheets=external_stylesheets)
    
    # Step 2: Register simulation and renderer callbacks
    set_callbacks!(app,
        initialize_sim,
        simulate,
        renderers,
        get_interactive_components(control_panel))

    # Step 3: Run the dashboard server
    run_server(app, host, port)

    return nothing
end
end # module ControlDashboard
