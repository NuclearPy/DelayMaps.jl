"""
    DelayMap{T<:Function}

Compactly supported linear operator with effective domain and codomain.

Fields:
- vector_field :: T

Constructors:
- `DelayMap(::Function)`

# Examples

```jldoctest
julia> DelayMap(sin)
DelayMap{typeof(sin)}(sin)
```
"""
struct DelayMap{T <: Function}
    vector_field :: T
end