# --- Example Usage ---
# The following code demonstrates how to use the module with a simple sine wave simulation.
using ControlDashboard, ControlDashboard.ControlPanel
using DataFrames
using Dash
using PlotlyJS

function run_sin_simulation((a, f, p, dt))
    @info "Running simulation"
    state = Dict(
        "amplitude" => a,
        "frequency" => f,
        "phase" => p,
        "duration" => 10,
        "dt" => dt,
    )
    t = 0:dt:10
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

# Create the app by passing our custom simulation functions to the dashboard wrapper
sinusoid_control_panel = make_panel(
    sin_wave_interfaces;
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

sinusoid_app =
    initialize_dashboard("Sinusoid Tuner"; control_panel = sinusoid_control_panel)
set_callbacks!(
    sinusoid_app,
    run_sin_simulation,
    Dict("main_view" => make_timeseries_figure),
    [
        ("amplitude", "value"),
        ("frequency", "value"),
        ("phase", "value"),
        ("dt", "value"),
    ],
)

# --- Main execution block ---
function main()
    # Run the server
    # You can access the dashboard at http://127.0.0.1:8050
    run_server(sinusoid_app, "0.0.0.0", 8050)
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
