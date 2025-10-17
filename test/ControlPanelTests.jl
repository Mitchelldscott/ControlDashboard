using Test, Dash, ControlDashboard

@testset "ControlPanel tests" begin
    # --- Test 1: input & slider with defaults ---
    @testset "input, slider and dropdown with defaults" begin
        comps = [
            Dict("component" => "input", "label" => "Initial Roll (deg)", "id" => "roll"),
            Dict("component" => "slider", "label" => "Pitch Angle", "id" => "pitch"),
            Dict(
                "component" => "dropdown",
                "label" => "Flight Mode",
                "id" => "mode",
                "options" => [
                    Dict("label" => "Stabilize", "value" => "stabilize"),
                    Dict("label" => "Acro", "value" => "acro"),
                    Dict("label" => "Loiter", "value" => "loiter"),
                ],
                "value" => "stabilize",
            ),
        ]
        out = make_panel(comps)
        @test out isa Component
        @test length(out.children) == 3

        # Test Input
        @test out.children[1].children[2].id == "roll"
        # Test Slider
        @test out.children[2].children[2].id == "pitch"
        # Test Dropdown
        dropdown_component = out.children[3].children[2]
        @test dropdown_component.id == "mode"
        @test dropdown_component.value == "stabilize"
        @test length(dropdown_component.options) == 3
        @test dropdown_component.options[1]["value"] == "stabilize"
    end

    # --- Test 2: input, slider and dropdown array ---
    @testset "input, slider and dropdown array" begin
        comps = [
            Dict(
                "component" => "input",
                "label" => "Roll",
                "id" => "roll",
            ),
            Dict(
                "component" => "input",
                "label" => "Pitch",
                "id" => "pitch",
            ),
            Dict(
                "component" => "slider",
                "label" => "Yaw",
                "id" => "yaw",
                "min" => -180,
                "max" => 180,
                "step" => 5,
            ),
            Dict(
                "component" => "dropdown",
                "label" => "Throttle Mode",
                "id" => "throttle_mode",
                "options" => [
                    Dict("label" => "Auto", "value" => "auto"),
                    Dict("label" => "Manual", "value" => "manual"),
                ],
                "value" => "auto",
            ),
        ]
        # a single container component
        grid = make_panel(comps)

        # The result 'grid' is an Array containing two sub Arrays.
        @test grid isa Component
        @test grid.id == "ControlPanel"
        @test length(grid.children) == 4

        # Now, test the children of the container.
        roll_input_container = grid.children[1]
        @test roll_input_container isa Component

        pitch_input_container = grid.children[2]
        @test pitch_input_container isa Component

        yaw_slider_container = grid.children[3]
        @test yaw_slider_container isa Component

        throttle_dropdown_container = grid.children[4]
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

        @test panel.id == "ControlPanel"
        @test length(panel.children) == 4

        component_1 = panel.children[1].children[2]
        @test component_1.id == "a"
        @test component_1.type == "number"

        component_2 = panel.children[2].children[2]
        @test component_2.id == "b"

        component_3 = panel.children[3].children[2]
        @test component_3.id == "c"
        @test component_3.type == "text"

        component_4 = panel.children[4].children[2]
        @test component_4.id == "d"
        @test component_4.data ==
              [Dict("d_1" => 1), Dict("d_1" => 2), Dict("d_1" => 3)]

        # Test custom shape
        panel2 = make_control_panel(params)
        @test length(panel2.children) == 4
        @test panel2.children[1].children[2].id == "a"
        @test panel2.children[2].children[2].id == "b"
        @test panel2.children[3].children[2].id == "c"
        @test panel2.children[4].children[2].id == "d"

        # Test fewer fields than shape
        small = SmallStruct(10, false)
        panel3 = make_control_panel(small)
        @test length(panel3.children) == 2
        @test panel3.children[1].children[2].id == "x"
        @test panel3.children[2].children[2].id == "y"

        # Test empty shape defaults
        panel4 = make_control_panel(small)
        @test length(panel4.children) == 2
    end
end

#--- get_interactive_components Test Suite ---

@testset "get_interactive_components Tests" begin
    @testset "Basic Cases" begin
        @testset "Handles empty panel" begin
            @test get_interactive_components(make_panel(Dict[])) == []
        end

        @testset "Extracts a single ID" begin
            config = [Dict("component" => "input", "label" => "Test", "id" => "test-id-1")]
            panel = make_panel(config)
            @test get_interactive_components(panel) == [("test-id-1", "value")]
        end

        @testset "Extracts multiple unique IDs" begin
            config = [
                Dict(
                    "component" => "input",
                    "label" => "A",
                    "id" => "id-a",
                ),
                Dict(
                    "component" => "input",
                    "label" => "B",
                    "id" => "id-b",
                ),
                Dict(
                    "component" => "datatable",
                    "label" => "C",
                    "id" => "id-c",
                ),
            ]
            panel = make_panel(config)
            expected_ids = [("id-a", "value"), ("id-b", "value"), ("id-c", "data")]
            results = get_interactive_components(panel)
            @test results == expected_ids
        end
    end

    @testset "Advanced Cases" begin
        @testset "Handles components with no ID" begin
            config = [
                Dict("component" => "slider", "label" => "A", "id" => "id-a"),
                Dict("component" => "checkbox", "label" => "B"), # No ID, wont register
                Dict("component" => "dropdown", "label" => "C", "id" => "id-c"),
            ]
            panel = make_panel(config)
            results = get_interactive_components(panel)
            @test results == [("id-a", "value"), ("id-c", "value")]
        end

        @testset "Returns only unique IDs when duplicates are present" begin
            config = [
                Dict("component" => "slider", "label" => "A", "id" => "id-a"),
                Dict("component" => "checkbox", "label" => "B", "id" => "id-b"),
                Dict("component" => "datatable", "label" => "C", "id" => "id-a"), # Duplicate ID
            ]
            panel = make_panel(config)
            results = get_interactive_components(panel)
            @test results == [("id-a", "value"), ("id-b", "value")]
        end

        @testset "Works on complex, nested panel structure" begin
            quadcopter_config = [
                Dict(
                    "component" => "input",
                    "label" => "Duration",
                    "id" => "t_final",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Sample time",
                    "id" => "dt",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "Roll",
                    "id" => "roll",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "Pitch",
                    "id" => "pitch",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "Yaw",
                    "id" => "yaw",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "P",
                    "id" => "p",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "Q",
                    "id" => "q",
                ),
                Dict(
                    "component" => "slider",
                    "label" => "R",
                    "id" => "r",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Kp",
                    "id" => "Kp",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Ki",
                    "id" => "Ki",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Kd",
                    "id" => "Kd",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Arm length",
                    "id" => "L",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Ixx",
                    "id" => "Ixx",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Iyy",
                    "id" => "Iyy",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Izz",
                    "id" => "Izz",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Mass",
                    "id" => "m",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Thrust Coeff",
                    "id" => "Kf",
                ),
                Dict(
                    "component" => "input",
                    "label" => "Drag Coeff",
                    "id" => "Km",
                ),
            ]
            expected = [(c["id"], "value") for c in quadcopter_config]
            panel = make_panel(quadcopter_config)
            results = get_interactive_components(panel)

            # The function should find all 18 unique IDs
            @test length(results) == 18
            @test sort(results) == sort(expected)
        end
    end
end
