# iterates through transiant orbits to find attracting periodic orbits
function transients(F::DelayMap, x::Sequence, τ::Real, iter::Int)
    a = copy(x)
    for n in 1:iter
        a = F(a,τ)
    end
    return a
end

# returns orbits of the system as represented by the coefficients of the chosen basis
function coefSamples(F::DelayMap, x0::Sequence, τ::Real, iter::Int)
    evals = zeros(dimension(space(x0)), iter)
    x = copy(x0)
    evals[:,1] = x.coefficients
    for i in 2:iter
        x = F(x, τ)
        evals[:,i] = x.coefficients
    end

    return evals
end

# returns orbits of the system as represented by evaluations of the continuous functions in orbits
# Times should be a range or array of times
function evalSamples(F::DelayMap, x0::Sequence, τ::Real, Times, iter::Int)
    evals = zeros(length(Times), iter)
    x = copy(x0)
    evals[:, 1] = x.(Times)
    for i in 2:iter
        x = F(x, τ)
        evals[:, i] = x.(Times)
    end

    return evals
end