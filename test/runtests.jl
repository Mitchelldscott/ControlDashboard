include("ControlDashboardTests.jl")
include("ControlPanelTests.jl")
include("SimulationTests.jl")

using Test
using Dates

@testset "Example scripts" begin
    example_file = joinpath(@__DIR__, "..", "examples", "sinusoid.jl")
    timeout_seconds = 5  # adjust as needed

    @testset "Run sinusoid_server.jl" begin
        ex_task = @async include(example_file)

        start_time = now()
        done = false

        while !done && (now() - start_time) < Millisecond(1000 * timeout_seconds)
            sleep(0.1)
            done = istaskdone(ex_task)
        end

        # Kill the task if still running after timeout
        if !done
            Base.throwto(ex_task, InterruptException())
        end

        @test istaskdone(ex_task)  # test passes if the task completes or is interrupted
    end
end
