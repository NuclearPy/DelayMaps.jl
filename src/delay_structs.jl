"""
    DelayMap{T <: Function}

Abstract type for all delay maps.
"""
abstract type DelayMap{ T<: Function} end

# Delay maps are functions from C[-1,1] -> C[-1,1], these constants hardcode this fact into the code
const PHASE_SPACE_INTERVAL_START = -1.0
const PHASE_SPACE_INTERVAL_END = 1.0
const PHASE_SPACE_INTERVAL_LENGTH = PHASE_SPACE_INTERVAL_END - PHASE_SPACE_INTERVAL_START

"""
    ExplicitDelayMap{T<:Function}

Delay map associated to a DDE with one constant delay, for which delay maps may be written explicitly.

Fields:
- vector_field :: T

Constructors:
- `ExplicitDelayMap(::Function)`

# Examples

```jldoctest
julia> ExplicitDelayMap(sin)
ExplicitDelayMap{typeof(sin)}(sin)
```
"""
struct ExplicitDelayMap{T <: Function} <: DelayMap
    vector_field :: T
end

function (F::ExplicitDelayMap)(τ::Real, x::Sequence)
    ∫fx = integrate(F.vector_field(x))
    ∫fx = ∫fx - ∫fx(PHASE_SPACE_INTERVAL_START)
    return x(PHASE_SPACE_INTERVAL_END) + τ * ∫fx / PHASE_SPACE_INTERVAL_LENGTH
end