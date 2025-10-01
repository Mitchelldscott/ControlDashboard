# ControlDashboard.jl
# A generalized Julia module to create an interactive control simulation dashboard using Dash.jl.
module ControlDashboard

    # Import required libraries
    using Dash
    using DataFrames
    
    include("ControlPanel.jl")

    # Export the main dashboard function so it can be used by other modules.
    export initialize_dashboard, set_callbacks

    """
        initialize_dashboard(title; interfaces=make_interfaces(), graphs=make_graphs()) -> DashApp

    Factory function that initializes a Dash app with the provided interfaces and graphs.
    """
    function initialize_dashboard(title; 
        interfaces=sample_time_and_duration_sliders(), 
        graphs=[dcc_graph(id = "main_view")]
    )
        app = dash(external_stylesheets=["https://bootswatch.com/5/cyborg/bootstrap.min.css"])

        app.layout = html_div() do
            html_h1(title, style = Dict("textAlign" => "center")),
            html_div(interfaces),  # user-defined or default interfaces
            graphs...               # user-defined or default graphs
        end

        return app
    end

    """
        control_dashboard(state_factory, run_simulation, make_figure)

    Create a generic Dash app for interactive simulation.

    # Arguments
    - `state_factory::Function`: A function that takes no arguments and returns a fresh 
    initial state object (e.g., a Dict or a custom struct with initial conditions).
    - `run_simulation::Function`: A function that executes the simulation. It should accept
    the initial state, duration, and dt as arguments and return the simulation data.
    - `make_figure::Function`: A function that takes the simulation data and returns a
    figure to be displayed.
    """
    function set_callbacks(app, state_factory, run_simulation, figures, interfaces)
        # Define the callback function that links the sliders to the graph.
        # When a slider value changes, this function is triggered.
        for (figure_name, renderer) in figures
            callback!(
                app,
                Output(figure_name, "figure"),
                [Input(interface, "value") for interface in interfaces]
            ) do interfaces...
                df = run_simulation(state_factory(interfaces))
                return renderer(df)
            end
        end
        # Return the configured app object
        return app
    end
end # module ControlDashboards