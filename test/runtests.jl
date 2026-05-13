# runtests.jl — test suite for BifurcationForecasting.jl
#
# Test cases:
#   1. Analytical 1-DOF scalar system (supercritical and subcritical Hopf)
#      Ref: Riso, Cesnik & Epureanu (2021) J. Fluids Struct. 101, §3
#   2. 4-state aeroelastic typical section (full, Riso 2022)
#      Ref: Riso, Cesnik & Epureanu (2022) AIAA J. 60, 5401–5413

using Test
using BifurcationForecasting

@testset "BifurcationForecasting.jl" begin
    include("test_analytical_1dof.jl")
    include("test_aeroelastic_full.jl")
end
