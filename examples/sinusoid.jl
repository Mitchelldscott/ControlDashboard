# --- Example Usage ---
# The following code demonstrates how to use the module with a simple sine wave simulation.
include("../src/ControlDashboard.jl")

using ControlDashboard
using ControlDashboard.ControlPanel
using DataFrames
using Dash
using PlotlyJS

# PlotlyBase.default_layout_template[] = "plotly_dark"

function run_sin_simulation((a, f, p, dt))
    @info "Running simulation"
    Dict(
        "amplitude" => a,
        "frequency" => f,
        "phase" => p,
        "duration" => 10,
        "dt" => dt,
    )
    t = 0:state["dt"]:state["duration"]
    values = state["amplitude"] .* sin.(state["frequency"] .* t .+ state["phase"])
    return DataFrame(; time = collect(t), value = values)
end

# 3. Figure Generator: Takes a DataFrame and creates a PlotlyJS plot.
function make_timeseries_figure(data)
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
        "component" => "slider",
        "label" => "Amplitude",
        "id" => "amplitude",
        "min" => 1.0,
        "max" => 100.0,
        "step" => 1.0,
        "value" => 1.0,
    ),
    Dict(
        "component" => "slider",
        "label" => "Frequency",
        "id" => "frequency",
        "min" => 1e-3,
        "max" => 10.0,
        "step" => 1e-2,
        "value" => 1.0,
    ),
    Dict(
        "component" => "slider",
        "label" => "Phase",
        "id" => "phase",
        "min" => -5.0,
        "max" => 5.0,
        "step" => 1e-1,
        "value" => 0.0,
    ),
    Dict(
        "component" => "slider",
        "label" => "Sample time",
        "id" => "dt",
        "min" => 1e-3,
        "max" => 1.0,
        "step" => 1e-2,
        "value" => 0.25,
    ),
]

# --- Main execution block ---
function main()
    # Create the app by passing our custom simulation functions to the dashboard wrapper
    panel = make_panel(
        sin_wave_interfaces;
        shape = (2, 2),
        component_style = Dict(
            "width" => "100%",        # Slider fills the available column
            "margin" => "4px 0",      # Vertical spacing between label and slider
            "display" => "block",
        ),
        label_style = Dict(
            "width" => "100%",
            "margin-bottom" => "4px",
            "font-weight" => "bold",
            "display" => "block",
            "text-align" => "center", # Center the label text
        ),
        panel_style = Dict(
            "width" => "100%",
            "display" => "grid",
            "grid-template-columns" => "1fr 1fr", # Two columns
            "gap" => "16px",                        # Space between columns/rows
            "align-items" => "stretch",             # Make children stretch full height
            "padding" => "16px",
        ),
    )
    app = initialize_dashboard("Sinusoid Tuner"; control_panel = panel)
    set_callbacks!(
        app,
        run_sin_simulation,
        Dict("main_view" => make_timeseries_figure),
        [
            ("amplitude", "value"),
            ("frequency", "value"),
            ("phase", "value"),
            ("dt", "value"),
        ],
    )

    # Run the server
    # You can access the dashboard at http://127.0.0.1:8050
    run_server(app, "0.0.0.0", 8050)
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
