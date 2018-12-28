"""
# module DataStructure



# Examples

```jldoctest
julia>
```
"""
module DataStructure

export Cell

struct Cell
    lattice::AbstractMatrix
    positions::AbstractMatrix
    numbers::AbstractVector
end

end