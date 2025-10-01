module ControlPanel
    
    using Dash

    export sample_time_and_duration_sliders, make_panel, build_component

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

    # builder for each component
    function build_component(spec::Dict, row_style=Dict())

        base_row_style = Dict("display"=>"flex", "alignItems"=>"center", "gap"=>"10px", "marginBottom"=>"8px")
        base_label_style = Dict("width"=>"160px", "marginRight"=>"6px")
        base_input_style = Dict("width"=>"140px")
        
        ctype = get(spec, "component", "input")  # default to input
        label = get(spec, "label", "")
        compid = get(spec, "id", "")

        component = nothing
        if ctype == "input"
            component = dcc_input(
                id=compid,
                type=get(spec,"type","number"),
                value=get(spec,"value",0.0),
                step=get(spec,"step",0.1),
                debounce=get(spec,"debounce",true),
                inputMode=get(spec,"inputmode","numeric"),
                min=get(spec,"min",nothing),
                max=get(spec,"max",nothing),
                style=base_input_style
            )
        elseif ctype == "slider"
            marks = Dict(string(i)=>i for i in get(spec,"min",0.0):get(spec,"step",0.1):get(spec,"max",1.0))
            component = dcc_slider(
                id=compid,
                min=get(spec,"min",0.0),
                max=get(spec,"max",1.0),
                step=get(spec,"step",0.1),
                value=get(spec,"value",0.0),
                marks=get(spec,"marks", marks),
                tooltip=Dict("always_visible"=>true)
            )
        elseif ctype == "dropdown"
            component = dcc_dropdown(
                id=compid,
                options=get(spec,"options",[]),
                value=get(spec,"value",nothing),
                clearable=get(spec,"clearable",true),
                searchable=get(spec,"searchable",true),
                multi=get(spec,"multi",false),
                style=base_input_style
            )
        else
            error("Unknown component type: $ctype")
        end

        return html_div([
            html_label(label, style=base_label_style),
            component
        ], style=merge(base_row_style, row_style))
    end

    """
    make_panel(configs::Vector{Dict}; shape=(nothing, nothing), row_style=Dict(), col_style=Dict())

    Build a flexible Dash panel of UI components.

    Each `spec` is a `Dict` describing a component. Minimum keys:
    - `"component"` :: String, one of `"input"`, `"slider"`, `"dropdown"`, ...
    - `"label"`     :: String, label text
    - `"id"`        :: String, component id

    Optional keys (depending on component type):
    - `"type"`      :: For input: `"number"`, `"text"`, ...
    - `"inputmode"` :: For input: `"decimal"`, `"numeric"`, ...
    - `"value"`     :: Default value
    - `"min"`, `"max"` :: Bounds
    - `"step"`      :: Step size
    - `"options"`   :: For dropdown: array of Dicts with `label` and `value`
    - `"marks"`     :: For slider: Dict of tick marks
    - `"debounce"`  :: For input: Bool
    - `"position"`  :: (row, col) grid position if `shape` provided

    Keywords:
    - `shape=(nrows, ncols)` → arrange components in a grid
    - `row_style`, `col_style` → CSS dicts merged with defaults
    """
    function make_panel(configs::Vector{<:Dict}; shape=(1,length(configs)),
                        row_style=Dict(), col_style=Dict())

        # vertical stack if no shape
        if shape == (1,length(configs))
            return [build_component(spec, row_style) for spec in configs]
        else
            nrows, ncols = shape
            grid = [html_div([]) for _ in 1:(nrows*ncols)]
            for spec in configs
                pos = get(spec, "position", nothing)
                if pos !== nothing
                    i,j = pos
                    idx = (i-1)*ncols + j
                    grid[idx] = build_component(spec, row_style)
                else
                    empty_idx = findfirst(x -> isempty(x.children), grid)
                    if empty_idx !== nothing
                        grid[empty_idx] = build_component(spec, row_style)
                    end
                end
            end
            rows = []
            for i in 1:nrows
                start = (i-1)*ncols + 1
                stop  = i*ncols
                push!(rows,
                    html_div(
                        grid[start:stop], 
                        style=merge(Dict("display"=>"flex","gap"=>"20px"), col_style)
                    )
                )
            end
            return rows
        end
    end
end # module
