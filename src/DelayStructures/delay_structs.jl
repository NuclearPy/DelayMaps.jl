"""
    DelayMap{T <: Function, S <: Function}

Abstract type for all delay maps.
"""
abstract type DelayMap{T<: Function, S <: Function} end

# Delay maps are functions from C[-1,1] -> C[-1,1], these constants hardcode this fact into the code
const PHASE_SPACE_INTERVAL_START = -1.0
const PHASE_SPACE_INTERVAL_END = 1.0
const PHASE_SPACE_INTERVAL_LENGTH = PHASE_SPACE_INTERVAL_END - PHASE_SPACE_INTERVAL_START

"""
    ExplicitDelayMap{T<:Function, S <: Function}

Delay map associated to a DDE with one constant delay, for which delay maps may be written explicitly. This struct is callable.

Fields:
- vector_field :: T
- derivative :: Sequence

Constructors:
- `ExplicitDelayMap(::Function, ::Function)`

# Examples

```jldoctest
julia> ExplicitDelayMap(sin, cos)
ExplicitDelayMap{typeof(sin), typeof(cos)}(sin, cos)
```
"""
struct ExplicitDelayMap{T <: Function, S <: Function} <: DelayMap
    vector_field :: T
    derivative :: S
end

# Solves a DDE with initial condition x as in the method-of-steps, translating time to [-1,1].
function (F::ExplicitDelayMap)(x::Sequence, τ::Real)
    ∫Fx = integrate(F.vector_field(x)); # ∫Fx = ∫Fx - ∫Fx(PHASE_SPACE_INTERVAL_START) <- include if not starting at -1
    return project(x(PHASE_SPACE_INTERVAL_END) + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x))
end

# Jacobian matrix of the delay map projected onto a finite basis
function Jacobian(F::ExplicitDelayMap, x::Sequence, τ::Real)
    dFx = F.derivative(x); dFx = project(Multiplication(dFx), space(x), image(Multiplication(dFx), space(x)))
    # E = project(Evaluation(PHASE_SPACE_INTERVAL_START), space(x), image(Multiplication(dFx), space(x)))
    ∫ = integralMatrix(codomain(dFx)); # ∫ = ∫ - E * ∫ <- include if not starting at -1
    E = project(Evaluation(1), space(x), codomain(∫))
    return E + τ*∫*dFx
end