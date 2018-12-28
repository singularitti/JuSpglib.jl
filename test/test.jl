module JuSpglibTests

using Test

using JuSpglib

@testset "Get symmetry of a simple lattices" begin
    cell = Cell([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0], transpose([0 0 0; 0.25 0.25 0.25]), [:Si, :Si])
    symops = get_symmetry(cell)
    @test length(symops) == 48
    symbol = get_international(cell)
    @test symbol == "Fd-3m"
    symbol = get_schoenflies(cell)
    @test symbol == "Oh^7"
end

@testset "Zincblende" begin
    cell = Cell([0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0], transpose([0 0 0; 0.25 0.25 0.25]), [:Zn, :S])
    symops = get_symmetry(cell)
    @test length(symops) == 24
    symbol = get_international(cell)
    @test symbol == "F-43m"
    symbol = get_schoenflies(cell)
    @test symbol == "Td^2"
end

end
