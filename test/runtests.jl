include("ControlDashboardTests.jl")
include("ControlPanelTests.jl")
include("SimulationTests.jl")

## Dry run examples
@testset "Sinusoid Example tests" begin
    include("../examples/sinusoid.jl")
    @testset "Panel Creation" begin
        panel = nothing
        @test_nowarn panel = make_panel(
            sin_wave_interfaces;
            component_style = Dict("width" => "100%"),
            label_style = Dict("font-weight" => "bold"),
            panel_style = Dict("display" => "grid"),
        )
        @test panel isa Component
        @test length(panel.children) == 4
    end
end

@testset "Quadcopter Attitude Stabilizer test" begin
    include("../examples/quadcopter_attitude_stabilization.jl")
    @testset "Panel Creation" begin
        panel = nothing
        @test_nowarn panel = quadcopter_interfaces()
        @test panel isa Component
        @test length(panel.children) == 18
    end
end
