module JuSpglibTests

using Test

using JuSpglib

@testset "Get symmetry of a simple lattices" begin
    cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
    positions = transpose([0 0 0; 0.25 0.25 0.25])
    types = [:Si, :Si]
    symops = get_symmetry(cell, positions, types)
    @test length(symops) == 48
    symbol = get_international(cell, positions, types)
    @test symbol == "Fd-3m"
    symbol = get_schoenflies(cell, positions, types)
    @test symbol == "Oh^7"
end

@testset "Zincblende" begin
    cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
    positions = transpose([0 0 0; 0.25 0.25 0.25])
    types = [:Zn, :S]
    symops = get_symmetry(cell, positions, types)
    @test length(symops) == 24
    symbol = get_international(cell, positions, types)
    @test symbol == "F-43m"
    symbol = get_schoenflies(cell, positions, types)
    @test symbol == "Td^2"
end

end
