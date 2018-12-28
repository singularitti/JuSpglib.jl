"""
# module DataStructure



# Examples

```jldoctest
julia>
```
"""
module DataStructure

export Cell

struct Cell{T <: AbstractFloat}
    lattice::AbstractMatrix{T}
    positions::AbstractMatrix{T}
    numbers::AbstractVector{T}
end

Cell(lattice::AbstractMatrix{T}, positions::AbstractMatrix{T}, numbers::AbstractVector{T}) where {T} = Cell{T}(lattice, positions, numbers)

end