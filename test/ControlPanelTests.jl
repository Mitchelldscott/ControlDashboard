using Test
using Dash
using ControlDashboard.ControlPanel

@testset "ControlPanel tests" begin
    # --- Test 1: input & slider with defaults ---
    @testset "input, slider and dropdown with defaults" begin
        comps = [
            Dict("component"=>"input", "label"=>"Initial Roll (deg)", "id"=>"roll"),
            Dict("component"=>"slider", "label"=>"Pitch Angle", "id"=>"pitch"),
            Dict(
                "component"=>"dropdown",
                "label"=>"Flight Mode",
                "id"=>"mode",
                "options"=>[
                    Dict("label"=>"Stabilize", "value"=>"stabilize"),
                    Dict("label"=>"Acro", "value"=>"acro"),
                    Dict("label"=>"Loiter", "value"=>"loiter"),
                ],
                "value"=>"stabilize",
            ),
        ]
        out = make_panel(comps)
        @test length(out) == 3
        # Test Input
        @test out[1].children[2].id == "roll"

        # Test Slider
        @test out[2].children[2].id == "pitch"

        # Test Dropdown
        dropdown_component = out[3].children[2]
        @test dropdown_component.id == "mode"
        @test dropdown_component.value == "stabilize"
        @test length(dropdown_component.options) == 3 # This should now pass
        @test dropdown_component.options[1]["value"] == "stabilize"
    end

    # --- Test 2: input, slider and dropdown array ---
    @testset "input, slider and dropdown array" begin
        comps = [
            Dict("component"=>"input", "label"=>"Roll", "id"=>"roll", "position"=>(1, 1)),
            Dict("component"=>"input", "label"=>"Pitch", "id"=>"pitch", "position"=>(1, 2)),
            Dict(
                "component"=>"slider",
                "label"=>"Yaw",
                "id"=>"yaw",
                "min"=>-180,
                "max"=>180,
                "step"=>5,
                "position"=>(2, 1),
            ),
            Dict(
                "component"=>"dropdown",
                "label"=>"Throttle Mode",
                "id"=>"throttle_mode",
                "options"=>[
                    Dict("label"=>"Auto", "value"=>"auto"),
                    Dict("label"=>"Manual", "value"=>"manual"),
                ],
                "value"=>"auto",
                "position"=>(2, 2),
            ),
        ]
        # a single container component
        grid = make_panel(comps; shape = (2, 2))

        # The result 'grid' is an Array containing two sub Arrays.
        @test grid isa Vector
        @test length(grid) == 2

        # Extract the rows from the array.
        row1 = grid[1]
        row2 = grid[2]
        @test row1 isa Component
        @test row2 isa Component

        # Now, test the children of the container.
        @test length(row1.children) == 2
        @test length(row2.children) == 2

        # Access components from the flat .children array
        roll_input_container = row1.children[1]
        pitch_input_container = row1.children[2]
        yaw_slider_container = row2.children[1]
        throttle_dropdown_container = row2.children[2]

        @test roll_input_container isa Component
        @test pitch_input_container isa Component
        @test yaw_slider_container isa Component
        @test throttle_dropdown_container isa Component

        # Test component IDs
        @test roll_input_container.children[2].id == "roll"
        @test pitch_input_container.children[2].id == "pitch"
        @test yaw_slider_container.children[2].id == "yaw"
        @test throttle_dropdown_container.children[2].id == "throttle_mode"

        # Test slider properties
        yaw_slider = yaw_slider_container.children[2]
        # @info yaw_slider
        @test yaw_slider.min == -180
        @test yaw_slider.max == 180
        @test yaw_slider.step == 5
    end

    # Test structs
    struct TestParams
        a::Int
        b::Bool
        c::String
        d::Vector{Int}
    end

    struct SmallStruct
        x::Int
        y::Bool
    end

    @testset "Control panel from struct tests" begin
        params = TestParams(42, true, "hello", [1, 2, 3])

        panel = make_control_panel(params)

        @test length(panel) == 4
        @test length(panel[1].children) == 2
        @test panel[1].children[2].id == "a"
        @test panel[1].children[2].type == "number"
        @test panel[2].children[2].id == "b"
        @test panel[3].children[2].id == "c"
        @test panel[3].children[2].type == "text"
        @test panel[4].children[2].id == "d"
        @test panel[4].children[2].options == [
            Dict("label" => "1", "value" => "1"),
            Dict("label" => "2", "value" => "2"),
            Dict("label" => "3", "value" => "3"),
        ]

        # Test custom shape
        panel2 = make_control_panel(params; shape = (2, 2))
        @test length(panel2) == 2
        @test length(panel2[1].children) == 2
        @test length(panel2[1].children[1].children) == 2
        @test panel2[1].children[1].children[2].id == "a"
        @test panel2[1].children[2].children[2].id == "b"
        @test panel2[2].children[1].children[2].id == "c"
        @test panel2[2].children[2].children[2].id == "d"

        # Test fewer fields than shape
        small = SmallStruct(10, false)
        panel3 = make_control_panel(small; shape = (2, 1))
        @test length(panel3) == 2
        @test panel3[1].children[1].children[2].id == "x"
        @test panel3[2].children[1].children[2].id == "y"

        # Test empty shape defaults
        panel4 = make_control_panel(small)
        @test length(panel4) == 2
        @test panel4[1].children[2].id == "x"
        @test panel4[2].children[2].id == "y"
    end
end
