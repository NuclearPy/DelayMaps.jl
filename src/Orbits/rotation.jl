# Averaging weight for computing rotation number given orbit data
function birkhoffWeight(t::Real)
    # Birkhoff weight function
    # t = 0 or t = 1 => w(t) = 0
    # 0 < t < 1 => w(t) = exp(1 / ((t - 1) * t))
    return (t == 0 || t == 1 ? 0 : exp(1 / ((t - 1) * t)))  # Birkhoff weight
end

function rotationNumberByEval(F::DelayMap, x0, τ, iter, center; freq = 2*π)
    A = 0
    sum = 0
    x = copy(x0)
#    firstInd = firstindex(x)
    for n in 0:(iter-1)
        y = x
        θ1 = mod2pi(atan(y(1) - center[2], y(-1)-center[1])) / freq
        x = F(x, τ)
        y = x
        θ2 = mod2pi(atan(y(1)-center[2], y(-1)-center[1]))/ freq

        wn = w(n/iter)
        A += wn
        pi2_Over_f = 2*π/freq
        sum += (θ2 < θ1 ? wn * (θ2 - θ1 + pi2_Over_f) : wn * (θ2 - θ1))
    end

    return sum/A
end

function rotationNumberByCoef(F::DelayMap, x0, τ, iter, center; freq = 2*π)
    A = 0
    sum = 0
    x = copy(x0)
#    firstInd = firstindex(x)
    for n in 0:(iter-1)
        y = x - center
        θ1 = mod2pi(atan(y[1], y[0])) / freq
        x = F(x, τ)
        y = x - center
        θ2 = mod2pi(atan(y[1], y[0]))/ freq

        wn = w(n/iter)
        A += wn
        pi2_Over_f = 2*π/freq
        sum += (θ2 < θ1 ? wn * (θ2 - θ1 + pi2_Over_f) : wn * (θ2 - θ1))
    end

    return sum/A
end