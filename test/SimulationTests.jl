include("../src/Simulation.jl")
include("../src/ControlDashboard.jl")
using Test, Dash, ControlDashboard, DifferentialEquations, StaticArrays, DataFrames

@testset "Simulation tests" begin
    @testset "ode_simulation tests" begin
        # Define a simple scalar ODE: du/dt = -u, solution u(t) = exp(-t)
        function simple_dynamics!(du, u, p, t)
            du[1] = -u[1]
        end

        # Run simulation
        df = rk4_simulation(
            simple_dynamics!,
            [1.0];
            t_final = 1.0,
            dt = 0.1,
            state_names = ["x1"],
        )

        # Check type
        @test df isa DataFrame

        # Check required columns
        @test names(df) == ["time", "x1"]

        # Time column goes from 0 → t_final with step dt
        @test df.time[1] ≈ 0.0
        @test df.time[end] ≈ 1.0

        # Analytical solution at t=1: exp(-1) ≈ 0.3679
        @test df.x1[end] ≈ exp(-1) atol = 1e-3
    end
end
