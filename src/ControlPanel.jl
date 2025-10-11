module ControlPanel

using Dash

export make_control_panel, make_panel, build_component, 
    get_interactive_components, sample_time_and_duration_sliders

"""
    build_component(config::Dict; component_style::Dict=Dict())

Create a standard Dash UI component (`input`, `slider`, or `dropdown`) with an associated label.  
This internal factory function standardizes how configuration dictionaries map to Dash components, simplifying UI assembly.

# Arguments
- `config::Dict`: **Required.** Defines the component type and its properties (see Configuration Keys below).
- `component_style::Dict`: **Optional.** Styles the outer `html_div` wrapper that contains the label and component. Defaults to `Dict()`.

# Returns
- `html_div`:
    - `html_label`: Label component centered above the interactive element.  
    - `element`: Interactive Dash HTML component created from the `config` dictionary.

# Configuration Keys
The `config` dictionary must contain the following foundational keys:

| Key           | Description |
|---------------|----------|
| `"component"` | Defines the component type: `"input"`, `"slider"`, or `"dropdown"`. |
| `"label"`     | The text displayed above the component. |
| `"id"`        | Unique component ID, required for Dash callbacks. |

---

## Component-Specific Keys

### **Input**
| Key | Default | Description |
|------|----------|-------------|
| `"type"`  | `"number"` | The HTML input type (`"text"`, `"number"`, etc.). |
| `"value"` | `0.0` | Initial input value. |
| `"step"`  | `1e-5` | Step increment for numeric inputs. |
| `"min"`   | `nothing` | Minimum allowable value. |
| `"max"`   | `nothing` | Maximum allowable value. |

### **Slider**
| Key | Default | Description |
|------|----------|-------------|
| `"min"`   | `0.0` | Minimum slider value. |
| `"max"`   | `1.0` | Maximum slider value. (Error if `min > max`.) |
| `"step"`  | `0.1` | Step increment. |
| `"marks"` | *(Calculated)* | Dict of marks; defaults to 5 evenly spaced values. |

### **Dropdown**
| Key | Default | Description |
|------|----------|-------------|
| `"options"` | `[]` | List of option dicts with `"label"` and `"value"` keys. |
| `"value"`   | `""` | Initially selected value(s). |
| `"multi"`   | `false` | Allow multiple selections if `true`. |

---

# Example

```julia
slider_config = Dict(
    "component" => "slider",
    "id" => "weight-slider",
    "label" => "Select Weight (kg)",
    "min" => 50.0,
    "max" => 100.0,
    "step" => 0.5,
    "value" => 75.0
)
```
"""
function build_component(config::Dict, component_style = Dict())

    ctype = get(config, "component", "input")
    label = get(config, "label", "")
    compid = config["id"]

    # Dispatch based on the component type string, using Val(Symbol(ctype))
    component = build_component_ui(Val(Symbol(ctype)), config, compid)

    # Common wrapper logic
    return html_div(
        [
            html_label(
                label,
                style = Dict(
                    "display" => "block",       # ensures label sits above
                    "text-align" => "center",   # center label over component
                    "margin-bottom" => "6px",
                    "font-weight" => "bold",
                ),
            ),
            component,
        ],
        style = component_style,
    )
end

# --- Specific Component Builders using Multiple Dispatch ---
# Note: dcc_input, dcc_slider, etc., and html_div, html_label are assumed to be defined by Dash

# Fallback method for unknown component types
function build_component_ui(::Val{S}, config::Dict, compid::AbstractString) where S
    error("Unknown component type: $S")
end

# 1. Input Component (Dispatches on ::Val{:input})
function build_component_ui(::Val{:input}, config::Dict, compid::AbstractString)
    return dcc_input(
        id = compid,
        type = get(config, "type", "number"),
        value = get(config, "value", 0.0),
        step = get(config, "step", 1e-7),
        debounce = get(config, "debounce", true),
        inputMode = get(config, "inputmode", "numeric"),
        min = get(config, "min", nothing),
        max = get(config, "max", nothing),
        style = Dict("width"=>"140px"),
    )
end

# 2. Slider Component (Dispatches on ::Val{:slider})
function build_component_ui(::Val{:slider}, config::Dict, compid::AbstractString)
    min_val = get(config, "min", 0.0)
    max_val = get(config, "max", 1.0)
    range = max_val - min_val

    if range < 0
        error("Slider Range must be greater than zero")
    end

    # Calculate default marks if none are provided
    marks = get(config, "marks", Dict(
        string(i)=>i for i = min_val:(range/5):max_val
    ))

    return dcc_slider(
        id = compid,
        min = min_val,
        max = max_val,
        step = get(config, "step", 0.1),
        value = get(config, "value", min_val),
        marks = marks,
        tooltip = Dict("always_visible"=>true),
    )
end

# 3. Dropdown Component (Dispatches on ::Val{:dropdown})
function build_component_ui(::Val{:dropdown}, config::Dict, compid::AbstractString)
    return dcc_dropdown(
        id = compid,
        options = get(config, "options", []),
        value = get(config, "value", ""),
        clearable = get(config, "clearable", true),
        searchable = get(config, "searchable", true),
        multi = get(config, "multi", false),
        style = Dict("width"=>"160px"),
    )
end

# 4. Checkbox Component (Dispatches on ::Val{:checkbox})
function build_component_ui(::Val{:checkbox}, config::Dict, compid::AbstractString)
    return dcc_checklist(
        id = compid,
        options = get(config, "options", []),
        value = get(config, "value", ""),
        style = Dict("width"=>"160px"),
    )
end

# 5. DataTable Component (Dispatches on ::Val{:datatable})
function build_component_ui(val::Val{:datatable}, config::Dict, compid::AbstractString)
    # Assumes DataTable is the Julia Dash component function name for dash_table.DataTable
    return dash_datatable(
        id = compid,
        columns = get(config, "columns", []),
        data = get(config, "data", val),
        editable = get(config, "editable", false),
        filter_action = get(config, "filter_action", "none"),
        sort_action = get(config, "sort_action", "none"),
        sort_mode = get(config, "sort_mode", "single"),
        column_selectable = get(config, "column_selectable", "none"),
        row_selectable = get(config, "row_selectable", "none"),
        row_deletable = get(config, "row_deletable", false),
        selected_columns = get(config, "selected_columns", []),
        selected_rows = get(config, "selected_rows", []),
        page_action = get(config, "page_action", "none"),
        page_current = get(config, "page_current", 0),
        page_size = get(config, "page_size", 10),
    )
end

"""
    make_panel(configs::Vector{Dict}; shape=(nothing, nothing), component_style=Dict(), col_style=Dict())

Build a flexible Dash panel of UI components.

Each `config` is a `Dict` describing a component. Required keys:
- `"component"` :: String, one of `"input"`, `"slider"`, `"dropdown"`, ...
- `"label"`     :: String, label text
- `"id"`        :: String, component id

Keywords:
- `shape=(nrows, ncols)` → arrange components in a grid
- `component_style`, `panel_style` → CSS dicts that specify the style of panels and components
"""
function make_panel(
    configs::Vector{<:Dict};
    shape = (1, length(configs)),
    component_style = Dict(),
    panel_style = Dict(),
)
    if shape[1] == 1
        return [build_component(config, component_style) for config in configs if "id" in keys(config)]
    else
        nrows, ncols = shape
        contents = [html_div([]) for _ = 1:(nrows*ncols)]
        for config in configs
            pos = get(config, "position", nothing)
            if !("id" in keys(config)) continue end
            if pos !== nothing
                i, j = pos
                idx = (i-1)*ncols + j
                contents[idx] = build_component(config, component_style)
            else
                empty_idx = findfirst(x -> isempty(x.children), contents)
                if empty_idx !== nothing
                    contents[empty_idx] = build_component(config, component_style)
                end
            end
        end
        rows = []
        for i = 1:nrows
            start = (i-1)*ncols + 1
            stop = i*ncols
            push!(rows, html_div(contents[start:stop], style = panel_style))
        end
        return rows
    end
end

# The main function is now just a dispatcher
"""
    infer_field_type!(config::Dict, val::Any) -> Dict

Infer the type of a Dash component based on the Julia value `val`
and update a configuration dictionary suitable for that UI component.
"""
function infer_field_type!(config::Dict, val)
    set_component_config!(config, val)
end

# --- Method for Booleans ---
function set_component_config!(config::Dict, val::Bool)
    config["component"] = "checkbox"
    config["value"] = val
end

# --- Method for Numbers ---
function set_component_config!(config::Dict, val::Number)
    config["component"] = "input"
    config["type"] = "number"
    config["value"] = val
end

# --- Method for Strings and Symbols (combined with Union) ---
function set_component_config!(config::Dict, val::Union{AbstractString, Symbol})
    config["component"] = "input"
    config["type"] = "text"
    config["value"] = string(val) # string() works for both
end

# --- Method for Vectors ---
function set_component_config!(config::Dict, val::AbstractVector)
    config["component"] = "datatable"
    
    # Define the columns: a single editable column for the vector elements
    config["columns"] = [
        Dict(
            "name" => "Value",
            "id" => "Value",
            "deletable" => true,
            "selectable" => true,
            "editable" => true # CRITICAL: Allows user to edit data
        )
    ]

    # Format the data: Convert the vector into a Vector of Dicts (list of rows)
    # Each original vector element becomes a row under the "Value" column.
    config["data"] = [
        Dict("Value" => o) for o in val
    ]

    # Set up standard DataTable features
    config["editable"] = true # Enable overall table editing
    config["filter_action"] = "native"
    config["sort_action"] = "native"
    config["row_deletable"] = true
    config["page_action"] = "native"
    config["page_current"] = 0
    config["page_size"] = 10
    config["column_selectable"] = "multi"
    config["row_selectable"] = "multi"
end

# --- Method for Enums ---
# Note: This is improved to correctly configure a dcc_slider
function set_component_config!(config::Dict, val::Enum)
    instances_list = instances(typeof(val))
    config["component"] = "slider"
    config["min"] = 1
    config["max"] = length(instances_list)
    config["step"] = 1
    # `marks` for a slider maps a numeric value to a string label
    config["marks"] = Dict(i => string(e) for (i, e) in enumerate(instances_list))
    # The slider's value will be the numeric index of the enum instance
    current_index = findfirst(isequal(val), instances_list)
    config["value"] = current_index
end

"""
    make_control_panel(params_struct::T; shape=nothing, component_style=Dict(), panel_style=Dict()) where {T}

Construct a 2D Dash control panel from a struct instance.

The function iterates through the fields of `params_struct`, creating a UI
component (a `dcc_input` box) for each field.

# Arguments
- `params_struct`: An instance of a Julia struct.
- `shape`: An optional `(rows, cols)` tuple to arrange the components in a grid.
        If not provided, it defaults to a single column.
- `component_style`, `panel_style`: Optional CSS style dictionaries passed to `make_panel`.

# Returns
- A `Vector` of `html_div` components representing the panel.
"""
function make_control_panel(
    params_struct::T;
    shape = nothing,
    component_style = Dict(),
    panel_style = Dict(),
) where {T}
    configs = Vector{Dict}()
    names = fieldnames(T)
    num_fields = length(names)

    # Default to a single column layout if no shape is provided
    final_shape = (shape === nothing) ? (1, num_fields) : shape
    nrows, ncols = final_shape

    # Check if the requested shape can hold all fields
    if num_fields > nrows * ncols
        @error "Shape ($nrows, $ncols) is too small for $num_fields fields."
    end

    for (i, name) in enumerate(names)
        val = getfield(params_struct, name)

        # Assign grid position based on row-major order
        row_pos = div(i - 1, ncols) + 1
        col_pos = mod(i - 1, ncols) + 1

        # Create the configuration dictionary for the component
        config = Dict(
            "id" => String(name),
            "label" => titlecase(String(name)),
            "position" => (row_pos, col_pos),
        )
        infer_field_type!(config, val)
        push!(configs, config)
    end

    # Use the provided `make_panel` function to build the final layout
    return make_panel(
        configs;
        shape = final_shape,
        component_style = component_style,
        panel_style = panel_style,
    )
end

"""
    sample_time_and_duration_sliders(; component_style::Dict, panel_style::Dict) -> Vector{DashComponent}

Default set of interface components (sliders + labels).
Users can supply their own instead.
"""
function sample_time_and_duration_sliders(component_style = Dict(), panel_style = Dict())
    comps = [
        Dict("component"=>"input", "label"=>"Sample Time", "id"=>"dt"),
        Dict("component"=>"input", "label"=>"Duration", "id"=>"t"),
    ]
    make_panel(comps; component_style = component_style, panel_style = panel_style)
end

function get_field_name(component_name)
    if component_name in ["dcc_input", "dcc_slider", "dcc_dropdown", "dcc_checklist"]
        return "value"
    elseif component_name == "dash_datatable"
        return "data"
    else
        return ""
    end
end

"""
    deduplicate_interfaces(interfaces)

Remove duplicate values from a list of tuples, only checks the first element.
"""
function deduplicate_interfaces(interfaces)
    seen = Set{String}()
    deduped = Tuple{String,String}[]
    for (id, field) in interfaces
        if !(id ∈ seen)
            push!(deduped, (id, field))
            push!(seen, id)
        end
    end
    return deduped
end

"""
    get_interactive_components(panel::Vector)

Return a vector of `(id, field)` pairs for all interactive components in the control panel.

This function:
- Recursively traverses nested `html_div` or `Component` elements.
- Identifies components with an `:id` property.
- Determines the correct field to read using `get_field_name(::Type)` dispatch.

Useful for building Dash callback interfaces automatically.
"""
function get_interactive_components(panel::Vector)
    interfaces = Tuple{String, String}[]
    for element in panel
        if element isa Component
            component_name = lowercase(String(getfield(element, 1)))
            # Collect (id, field) pairs for interactive components
            if hasproperty(element, :id) && !isnothing(element.id) && !isempty(element.id)
                id = String(element.id)
                field = get_field_name(component_name)
                if length(field) > 0
                    push!(interfaces, (id, field))
                end
            end
            # Recurse through children if present
            if hasproperty(element, :children) && !isempty(element.children)
                children = element.children isa AbstractVector ? collect(Iterators.flatten([element.children])) : []
                child_info = get_interactive_components(children)
                if length(child_info) > 0 
                    append!(interfaces, child_info)
                end
            end
        end
    end
    
    return deduplicate_interfaces(interfaces)
end

end # module
