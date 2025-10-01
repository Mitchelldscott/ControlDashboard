module ControlPanel
    
    export sample_time_and_duration_sliders

    """
    make_interfaces() -> Vector{DashComponent}

    Default set of interface components (sliders + labels).
    Users can supply their own instead.
    """
    function sample_time_and_duration_sliders()
        return [
            html_label("Simulation Duration"),
            dcc_slider(
                id = "duration",
                min = 1e-5,
                max = 100,
                step = 10,
                value = 50,
                marks = Dict([i => string(i) for i in 0:10:100])
            ),

            html_label("Sampling time (dt)"),
            dcc_slider(
                id = "dt",
                min = 1e-1,
                max = 10,
                step = 1e-1,
                value = 1,
                marks = Dict([i => string(i) for i in 0:1:10])
            ),
        ]
    end
end