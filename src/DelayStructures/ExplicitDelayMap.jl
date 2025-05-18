"""
    ExplicitDelayMap{T<:Function, S <: Function}

Delay map associated to a DDE with one constant delay, for which delay maps may be written explicitly. This struct is callable.
One dimensional DDE for now.

Fields:
- vector_field :: T
- derivative :: S

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

# --- Truncated delay map on Chebyshev basis ---
# Chebyshev forms the phase space of our system

# Solves a DDE with initial condition x as in the method-of-steps, translating time to [-1,1].
function (F::ExplicitDelayMap)(x::Sequence{Chebyshev, Vector{S}}, τ::Real) where S <: Real
    ∫Fx = integrate(F.vector_field(x)); # ∫Fx = ∫Fx - ∫Fx(PHASE_SPACE_INTERVAL_START) <- include if not starting at -1
    return project(x(PHASE_SPACE_INTERVAL_END) + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x))
end

# Jacobian matrix of the delay map projected onto a finite basis
function Jacobian(F :: ExplicitDelayMap, x :: Sequence{Chebyshev, Vector{S}}, τ :: Real) where S <: Real
    dFx = F.derivative(x); dFx = project(Multiplication(dFx), space(x), image(Multiplication(dFx), space(x)))
    # E = project(Evaluation(PHASE_SPACE_INTERVAL_START), space(x), image(Multiplication(dFx), space(x)))
    ∫ = integralMatrix(codomain(dFx)); # ∫ = ∫ - E * ∫ <- include if not starting at -1
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END), space(x), codomain(∫))
    return E + (τ / PHASE_SPACE_INTERVAL_LENGTH) * ∫ * dFx 
end

# --- Truncated delay map on Chebyshev ⊗ Fourier basis ---
# used for computing periodic orbits of DDE // quasi-periodic orbits of delay map

function (F::ExplicitDelayMap)(x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}}}, S}, τ::Real) where S <: AbstractArray
    ∫Fx = Integral(1,0) * F.vector_field(x)
    return project(x(PHASE_SPACE_INTERVAL_END, nothing) + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x))
end

function Jacobian(F::ExplicitDelayMap, x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}}}, S}, τ::Real) where S <: AbstractArray
    dfx = F.derivative(x)
    dfx = project(Multiplication(dfx), space(x), image(Multiplication(dfx), space(x)))
    ∫ = project(Integral(1,0), codomain(dfx), image(Integral(1,0), codomain(dfx)))
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END, nothing), space(x), codomain(∫))
    return project(E + τ*∫*dfx/PHASE_SPACE_INTERVAL_LENGTH, space(x), space(x))
end

# --- Truncated delay map on Chebyshev ⊗ Taylor basis ---
# used for computing manifolds at fixed points of delay map

# 1-dimensional manifolds
function (F::ExplicitDelayMap)(x::Sequence{TensorSpace{Tuple{Chebyshev, Taylor}}, S}, τ::Real) where S <: AbstractArray
    ∫Fx = Integral(1,0) * F.vector_field(x)
    return project(x(PHASE_SPACE_INTERVAL_END, nothing) + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x))
end

function Jacobian(F::ExplicitDelayMap, x::Sequence{TensorSpace{Tuple{Chebyshev, Taylor}}, S}, τ::Real) where S <: AbstractArray
    dfx = F.derivative(x)
    dfx = project(Multiplication(dfx), space(x), image(Multiplication(dfx), space(x)))
    ∫ = project(Integral(1,0), codomain(dfx), image(Integral(1,0), codomain(dfx)))
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END, nothing), space(x), codomain(∫))
    return project(E + τ*∫*dfx/PHASE_SPACE_INTERVAL_LENGTH, space(x), space(x))
end

# n-dimensional manifolds (n = Morse)
function (F::ExplicitDelayMap)(x::Sequence{TensorSpace{Tuple{Chebyshev, Vararg{Taylor, Morse}}}, S}, τ::Real) where S <: AbstractArray where Morse
    T = Sequence{TensorSpace{Tuple{Chebyshev, Vararg{Taylor, Morse}}}, S}
    ∫Fx = Integral(1, Tuple(0 for n in 1:Morse)...) * F.vector_field(x)
    Evalx = Evaluation(PHASE_SPACE_INTERVAL_END, Tuple(nothing for n in 1:Morse)...) * x
    return convert(T, project(Evalx + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x)))
end

function Jacobian(F::ExplicitDelayMap, x::Sequence{TensorSpace{Tuple{Chebyshev, Vararg{Taylor, Morse}}}, S}, τ::Real) where S <: AbstractArray where Morse
    T = LinearOperator{TensorSpace{Tuple{Chebyshev, Vararg{Taylor, Morse}}}, TensorSpace{Tuple{Chebyshev, Vararg{Taylor, Morse}}}, Matrix{eltype(S)}}
    dfx = F.derivative(x)
    dfx = project(Multiplication(dfx), space(x), image(Multiplication(dfx), space(x)))
    ∫ = project(Integral(1, Tuple(0 for n in 1:Morse)...), codomain(dfx), image(Integral(1, Tuple(0 for n in 1:Morse)...), codomain(dfx)))
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END, Tuple(nothing for n in 1:Morse)...), space(x), codomain(∫))
    return convert(T, project(E + τ*∫*dfx/PHASE_SPACE_INTERVAL_LENGTH, space(x), space(x)))
end

# --- Truncated delay map on Chebyshev ⊗ Fourier ⊗ Taylor basis ---
# used for computing manifolds at periodic orbits of DDE


function (F::ExplicitDelayMap)(x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Taylor}}, S}, τ::Real) where S <: AbstractArray
    ∫Fx = Integral(1,0,0) * F.vector_field(x)
    return project(x(PHASE_SPACE_INTERVAL_END, nothing, nothing) + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x))
end

function Jacobian(F::ExplicitDelayMap, x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Taylor}}, S}, τ::Real) where S <: AbstractArray
    dfx = F.derivative(x)
    dfx = project(Multiplication(dfx), space(x), image(Multiplication(dfx), space(x)))
    ∫ = project(Integral(1,0,0), codomain(dfx), image(Integral(1,0,0), codomain(dfx)))
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END, nothing, nothing), space(x), codomain(∫))
    return project(E + τ*∫*dfx/PHASE_SPACE_INTERVAL_LENGTH, space(x), space(x))
end

function (F::ExplicitDelayMap)(x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Vararg{Taylor, Morse}}}, S}, τ::Real) where S <: AbstractArray where Morse
    T = Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Vararg{Taylor, Morse}}}, S}
    ∫Fx = Integral(1, 0, Tuple(0 for n in 1:Morse)...) * F.vector_field(x)
    Evalx = Evaluation(PHASE_SPACE_INTERVAL_END, nothing, Tuple(nothing for n in 1:Morse)...) * x
    return convert(T, project(Evalx + τ * ∫Fx / PHASE_SPACE_INTERVAL_LENGTH, space(x)))
end

function Jacobian(F::ExplicitDelayMap, x::Sequence{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Vararg{Taylor, Morse}}}, S}, τ::Real) where S <: AbstractArray where Morse
    T = LinearOperator{TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Vararg{Taylor, Morse}}}, TensorSpace{Tuple{Chebyshev, Fourier{Float64}, Vararg{Taylor, Morse}}}, Matrix{eltype(S)}}
    dfx = F.derivative(x)
    dfx = project(Multiplication(dfx), space(x), image(Multiplication(dfx), space(x)))
    ∫ = project(Integral(1, 0, Tuple(0 for n in 1:Morse)...), codomain(dfx), image(Integral(1, 0, Tuple(0 for n in 1:Morse)...), codomain(dfx)))
    E = project(Evaluation(PHASE_SPACE_INTERVAL_END, nothing, Tuple(nothing for n in 1:Morse)...), space(x), codomain(∫))
    return convert(T, project(E + τ*∫*dfx/PHASE_SPACE_INTERVAL_LENGTH, space(x), space(x)))
end