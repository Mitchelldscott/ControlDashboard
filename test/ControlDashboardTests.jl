using Test, Dash, ControlDashboard, DifferentialEquations, StaticArrays, DataFrames

@testset "ControlDashboard tests" begin

    # --- 1. Setup Mock Environment ---
    # Define mock data structures that will be returned by our mock functions.
    MOCK_STATE = Dict("status" => "state_created")
    MOCK_DATAFRAME = (colA = [10, 20], colB = [30, 40]) # Using a NamedTuple as a mock DataFrame
    MOCK_FIGURE = Dict("data" => "mock_plot", "layout" => "mock_layout")

    # Use Refs to act as simple logs, tracking calls and arguments to the mock functions.
    state_factory_calls = Ref{Vector{Any}}([])
    run_simulation_calls = Ref{Vector{Any}}([])
    renderer_calls = Ref{Vector{Any}}([])

    # Define mock versions of the simulation and rendering functions.
    # Each function logs its input arguments before returning a predefined mock object.
    function mock_state_factory(inputs)
        push!(state_factory_calls[], inputs)
        return MOCK_STATE
    end

    function mock_run_simulation(state)
        push!(run_simulation_calls[], state)
        return MOCK_DATAFRAME
    end

    function mock_renderer(df)
        push!(renderer_calls[], df)
        return MOCK_FIGURE
    end

    # Helper function to clear the logs before each new testset.
    function reset_mocks()
        empty!(state_factory_calls[])
        empty!(run_simulation_calls[])
        empty!(renderer_calls[])
    end

    # --- 2. Define Test Cases ---

    @testset "Single figure and multiple interfaces" begin
        reset_mocks()

        # Define the inputs for this test case
        figures = Dict("main_graph" => mock_renderer)
        interfaces = [("duration_slider", "value"), ("dt_slider", "value")]

        # Execute the function to be tested
        set_callbacks!(dash(), mock_state_factory, mock_run_simulation, figures, interfaces)

        # --- 3. Assertions ---
        test_slider_values = (100, 0.5)
        df = mock_run_simulation(mock_state_factory(test_slider_values))
        result_figure = mock_renderer(df)

        # Verify the chain of function calls
        @test length(state_factory_calls[]) == 1
        @test first(state_factory_calls[]) == test_slider_values

        @test length(run_simulation_calls[]) == 1
        @test first(run_simulation_calls[]) == MOCK_STATE

        @test length(renderer_calls[]) == 1
        @test first(renderer_calls[]) == MOCK_DATAFRAME

        # Verify the final output
        @test result_figure == MOCK_FIGURE
    end

    @testset "Multiple figures setup" begin
        reset_mocks()

        app = dash()
        figures = Dict("figure_A" => mock_renderer, "figure_B" => mock_renderer)
        interfaces = [("sliderA", "value"), ("sliderB", "value")]

        set_callbacks!(app, mock_state_factory, mock_run_simulation, figures, interfaces)

        # Two callbacks should have been created
        @test length(app.callbacks) == 2
    end

    @testset "Edge case: No figures" begin
        reset_mocks()

        app = dash()
        figures = Dict() # Empty dictionary
        interfaces = [("slider", "value")]

        set_callbacks!(app, mock_state_factory, mock_run_simulation, figures, interfaces)

        # No callbacks should be registered if there are no figures
        @test isempty(app.callbacks)

        # None of the processing functions should have been called
        @test isempty(state_factory_calls[])
        @test isempty(run_simulation_calls[])
        @test isempty(renderer_calls[])
    end
end