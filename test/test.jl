module JuSpglibTests

using Test

using JuSpglib

@testset "Get symmetry of a simple lattices" begin
    cell = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]
    positions = transpose([0 0 0; 0.25 0.25 0.25])
    types = [:Si, :Si]
    symops = symmetry_operations(cell, positions, types)
    @test length(symops) == 48
    symbol = international_symbol(cell, positions, types)
    @test symbol == "Fd-3m"
    symbol = schoenflies_symbol(cell, positions, types)
    @test symbol == "Oh^7"
end

@testset "Zincblende" begin
    types = [:Zn, :S]
    symops = symmetry_operations(cell, positions, types)
    @test length(symops) == 24
    symbol = international_symbol(cell, positi®ons, types)
    @test symbol == "F-43m"
    symbol = schoenflies_symbol(cell, positions, types)
    @test symbol == "Td^2"
end

end
