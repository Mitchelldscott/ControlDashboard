# --- Example Usage ---
# The following code demonstrates how to use the module with a simple sine wave simulation.
using ControlDashboard
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
        "dt" => interfaces[4]
    )
end

function run_sin_simulation(state)
    @info "Running simulation"
    t = 0:state["dt"]:state["duration"]
    values = state["amplitude"] .* sin.(state["frequency"] .* t .+ state["phase"])
    return DataFrame(time = collect(t), value = values)
end

# 3. Figure Generator: Takes a DataFrame and creates a PlotlyJS plot.
function make_timeseries_figure(data::DataFrame)
    return Plot(
        scatter(x = data.time, y = data.value, mode = "lines", name="Signal"),
        Layout(
            title = "Sinusoid",
            xaxis_title = "Time (seconds)",
            yaxis_title = "Output",
            template="plotly_dark"
        )
    )
end

sin_wave_interfaces = [
    html_label("Amplitude"),
    dcc_slider(
        id = "amplitude", min = 1e-5, max = 10,
        step = 1, value = 1,
        marks = Dict([i => string(i) for i in 0:1:10])
    ),
    html_label("Frequency"),
    dcc_slider(
        id = "frequency", min = 1e-2, max = 10,
        step = 1e-2, value = 1,
        marks = Dict([i => string(i) for i in 0:1:10])
    ),
    html_label("Phase"),
    dcc_slider(
        id = "phase", min = -5, max = 5,
        step = 1e-1, value = 0,
        marks = Dict([i => string(i) for i in -5:1:5])
    ),
    html_label("Sample time"),
    dcc_slider(
        id = "dt", min = 1e-5, max = 10,
        step = 1e-1, value = 1,
        marks = Dict([i => string(i) for i in 0:1:10])
    ),
]


# --- Main execution block ---
function main()    
    println("Starting control dashboard...")

    # Create the app by passing our custom simulation functions to the dashboard wrapper
    app = initialize_dashboard("Sinusoid Tuner"; interfaces=sin_wave_interfaces)
    app = set_callbacks(
        app,
        initial_state,
        run_sin_simulation,
        Dict("main_view" => make_timeseries_figure),
        ["amplitude", "frequency", "phase", "dt"]
    )
    
    # Run the server
    # You can access the dashboard at http://127.0.0.1:8050
    run_server(app, "0.0.0.0", 8050)
end

# Run the main function
main()