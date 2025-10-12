# --- Example Usage ---
# The following code demonstrates how to use the module with a simple sine wave simulation.
using ControlDashboard
using ControlDashboard.ControlPanel
using DataFrames
using Dash
using PlotlyJS

# PlotlyBase.default_layout_template[] = "plotly_dark"

function initial_state(interfaces)
    return Dict(
        "amplitude" => interfaces[1],
        "frequency" => interfaces[2],
        "phase" => interfaces[3],
        "duration" => 10,
        "dt" => interfaces[4],
    )
end

function run_sin_simulation(state)
    @info "Running simulation"
    t = 0:state["dt"]:state["duration"]
    values = state["amplitude"] .* sin.(state["frequency"] .* t .+ state["phase"])
    return DataFrame(; time = collect(t), value = values)
end

# 3. Figure Generator: Takes a DataFrame and creates a PlotlyJS plot.
function make_timeseries_figure(data::DataFrame)
    @info "Rendering Plot"
    return Plot(
        scatter(; x = data.time, y = data.value, mode = "lines", name = "Signal"),
        Layout(;
            title = "Sinusoid",
            xaxis_title = "Time (seconds)",
            yaxis_title = "Output",
            template = "plotly_dark",
        ),
    )
end

sin_wave_interfaces = [
    Dict(
        "component"=>"slider",
        "label"=>"Amplitude",
        "id"=>"amplitude",
        "min"=>1.0,
        "max"=>10.0,
        "step"=>1.0,
        "value"=>1.0,
    ),
    Dict(
        "component"=>"slider",
        "label"=>"Frequency",
        "id"=>"frequency",
        "min"=>1e-2,
        "max"=>10.0,
        "step"=>1e-2,
        "value"=>1.0,
    ),
    Dict(
        "component"=>"slider",
        "label"=>"Phase",
        "id"=>"phase",
        "min"=>-5.0,
        "max"=>5.0,
        "step"=>1e-1,
        "value"=>0.0,
    ),
    Dict(
        "component"=>"slider",
        "label"=>"Sample time",
        "id"=>"dt",
        "min"=>1e-2,
        "max"=>1.0,
        "step"=>1e-2,
        "value"=>0.1,
    ),
]

# --- Main execution block ---
function main()
    # Create the app by passing our custom simulation functions to the dashboard wrapper
    panel = make_panel(
        sin_wave_interfaces;
        shape = (2, 2),
        panel_style = Dict(
            "width" => "100%",
            "display" => "flex",
            "flex-direction" => "row",
            "align-items" => "stretch",
            "justify-content" => "space-evenly",
        ),
        component_style = Dict(
            "flex" => "1",
            "width" => "100%",
            "display" => "flex",
            "flex-direction" => "column",
            "align-items" => "stretch",
        ),
    )
    app = initialize_dashboard("Sinusoid Tuner"; control_panel = panel)
    # set_callbacks!(
    #     app,
    #     initial_state,
    #     run_sin_simulation,
    #     Dict("main_view" => make_timeseries_figure),
    #     [("amplitude", "value"), ("frequency", "value"), ("phase", "value"), ("dt", "value")],
    # )

    # Run the server
    # You can access the dashboard at http://127.0.0.1:8050
    run_server(app, "0.0.0.0", 8050)
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
