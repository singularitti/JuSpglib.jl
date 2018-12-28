"""
# module CAPIs



# Examples

```jldoctest
julia>
```
"""
module CAPIs

using CoordinateTransformations
using DataStructures: counter

using JuSpglib.DataStructure: Cell

export get_symmetry, get_international, get_schoenflies

function getfields(obj, fields...)
    Tuple(getfield(obj, name) for name in fields)
end

function get_ccell(cell::Cell)::Cell
    lattice, positions, types = getfields(cell, :lattice, :positions, :numbers)
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    ctypes = convert(Vector{Cint}, [repeat([i], v) for (i, v) in (enumerate ∘ values ∘ counter)(types)] |> Iterators.flatten |> collect)
    return Cell(clattice, cpositions, ctypes)
end

cchars_to_string(s::Vector{Cchar}) = map(Char, s) |> join |> x -> split(x, "\0") |> first

function get_symmetry(cell::Cell; symprec::Real = 1e-8)
    size(cell.positions, 2) != length(cell.numbers) && throw(DimensionMismatch("The number of positions and atomic types do not match!"))
    size(cell.positions, 1) != 3 && error("Operations in 3D space is supported here!")

    maxsize = 52
    rotations = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    ccell = get_ccell(cell)
    clattice, cpositions, ctypes = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_symmetry, spglib), Cint,
        (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        rotations, translations, maxsize, clattice, cpositions, ctypes, length(ctypes), symprec)
    numops == 0 && error("Could not determine symmetries!")

    [AffineMap(transpose(rotations[:, :, i]), translations[:, i]) for i in 1:numops]
end

function get_international(cell::Cell; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    ccell = get_ccell(cell)
    clattice, cpositions, ctypes = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_international, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, ctypes, length(ctypes), symprec)
    numops == 0 && error("Could not determine the international symbol!")

    cchars_to_string(result)
end

function get_schoenflies(cell::Cell; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    ccell = get_ccell(cell)
    clattice, cpositions, ctypes = getfields(ccell, :lattice, :positions, :numbers)

    numops = ccall((:spg_get_schoenflies, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, ctypes, length(ctypes), symprec)
    numops == 0 && error("Could not determine the Schoenflies symbol!")

    cchars_to_string(result)
end

end