module ControlPanel
    
    using Dash

    export make_panel, build_component, sample_time_and_duration_sliders

    # builder for each component
    function build_component(config::Dict, row_style=Dict())

        base_row_style = Dict("display"=>"flex", "alignItems"=>"center", "gap"=>"10px", "marginBottom"=>"8px")
        base_label_style = Dict("width"=>"160px", "marginRight"=>"6px")
        base_input_style = Dict("width"=>"140px")
        
        ctype = get(config, "component", "input")  # default to input
        label = get(config, "label", "")
        compid = get(config, "id", "")

        component = nothing
        if ctype == "input"
            component = dcc_input(
                id=compid,
                type=get(config,"type","number"),
                value=get(config,"value",0.0),
                step=get(config,"step",1e-5),
                debounce=get(config,"debounce",true),
                inputMode=get(config,"inputmode","numeric"),
                min=get(config,"min",nothing),
                max=get(config,"max",nothing),
                style=base_input_style
            )
        elseif ctype == "slider"
            range = get(config,"max",1.0) - get(config,"min",0.0)
            if range < 0 error("Slider Range must be greater than zero") end

            marks = Dict(string(i)=>i for i in get(config,"min",0.0):range/5:get(config,"max",1.0))
            component = dcc_slider(
                id=compid,
                min=get(config,"min",0.0),
                max=get(config,"max",1.0),
                step=get(config,"step",0.1),
                value=get(config,"value",get(config,"min",0.0)), # Default must be valid
                marks=get(config,"marks",marks),
                tooltip=Dict("always_visible"=>true)
            )
        elseif ctype == "dropdown"
            component = dcc_dropdown(
                id=compid,
                options=get(config,"options",[]),
                value=get(config,"value",nothing),
                clearable=get(config,"clearable",true),
                searchable=get(config,"searchable",true),
                multi=get(config,"multi",false),
                style=base_input_style
            )
        else
            error("Unknown component type: $ctype")
        end

        return html_div([
            html_label(label, style=base_label_style),
            component
        ], style=base_row_style)
    end

    """
    make_panel(configs::Vector{Dict}; shape=(nothing, nothing), row_style=Dict(), col_style=Dict())

    Build a flexible Dash panel of UI components.

    Each `config` is a `Dict` describing a component. Minimum keys:
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
            return [build_component(config, row_style) for config in configs]
        else
            nrows, ncols = shape
            grid = [html_div([]) for _ in 1:(nrows*ncols)]
            for config in configs
                pos = get(config, "position", nothing)
                if pos !== nothing
                    i,j = pos
                    idx = (i-1)*ncols + j
                    grid[idx] = build_component(config, row_style)
                else
                    empty_idx = findfirst(x -> isempty(x.children), grid)
                    if empty_idx !== nothing
                        grid[empty_idx] = build_component(config, row_style)
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

    """
    make_interfaces() -> Vector{DashComponent}

    Default set of interface components (sliders + labels).
    Users can supply their own instead.
    """
    function sample_time_and_duration_sliders(row_style=Dict(), col_style=Dict())
        comps = [
            Dict("component"=>"input", "label"=>"Sample Time", "id"=>"dt"),
            Dict("component"=>"input", "label"=>"Duration", "id"=>"t"),
        ]
        make_panel(comps; row_style=row_style, col_style=col_style)
    end
end # module
