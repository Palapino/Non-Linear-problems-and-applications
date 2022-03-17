#=
Author : Mathis Antonetti
Date : 17/03/2022
Description :
This is a script to compare three analytical methods to tackle the vertical launching problem.
The non-friction case, the linear friction case, and the quadratic friction (viscosity of a newtonian fluid) case.
=#
using Pkg
using Plots

# General tools
function newton(f, fprime, z0, tol)
    u = z0
    n = 0
    while(abs(f(u)) > tol && n < 50)
        u = u - f(u)/fprime(u)
        n = n + 1
    end
    return u
end

function maxi(V::Vector{Float64})
    v = V[1]
    l = size(V)[1]
    for k=1:l
        v = max(v, V[k])
    end
    return v
end

# test for newton
begin
    function f_test(x)
        return exp(4*x) - (1 + x)
    end
    function f_test_prime(x)
        return 4*exp(4*x) - 1
    end
    @assert abs(f_test(newton(f_test, f_test_prime, 0.5, 1e-8))) <= 1e-8
end

## without friction :

# the time when the ball is at the highest position
function T_sf(v0, g)
    return v0/g
end

# the ball velocity
function v_sf(v0, g, t)
    return g*(T_sf(v0, g) - t)
end

# the ball velocity when it hit the floor
function vfall_sf(v0, g, h)
    return sqrt(2g*h + v0^2)
end

# the ball position
function x_sf(v0, g, h, t)
    return h + v0*t - 0.5*g*t^2
end

## linear friction

# the time when the ball is at the highest position
function T_ff(v0, g, beta)
    return 1/beta*log(1 + beta/g*v0)
end

# the ball velocity
function v_ff(v0, g, beta, t)
    return g/beta*(exp(beta*(T_ff(v0, g, beta) - t)) - 1)
end

# the ball velocity when it hit the floor
function vfall_ff(v0, g, h, beta)
    function _f(x)
        return log((1+beta/g*x)/(1+beta/g*v0)) - beta/g*(x - v0) - beta^2/g*h
    end
    function _fprime(x)
        return (1/(1+beta/g*x) - 1)*beta/g
    end
    return newton(_f, _fprime, v0/2.0, 1e-10)
end

# the ball position
function x_ff(v0, g, h, beta, t)
    return h - g/beta*t + 1/beta*(v0 + g/beta)*(1 - exp(-beta*t))
end

## quadratic friction

# the time when the ball is at the highest position
function T_vf(v0, g, alpha)
    return 1/sqrt(alpha*g)*atan(sqrt(alpha/g)*v0)
end

# the ball velocity
function v_vf(v0, g, alpha, t)
    Tfv = T_vf(v0, g, alpha)
    if(Tfv > t)
        return sqrt(g/alpha)*tan(sqrt(alpha*g)*(Tfv - t))
    else
        return sqrt(g/alpha)*tanh(sqrt(alpha*g)*(Tfv - t))
    end
end

# the ball position
function x_vf(v0, g, h, alpha, t)
    Tfv = T_vf(v0, g, alpha)
    if(Tfv > t)
        return 1/alpha*log(abs(cos(sqrt(alpha*g)*(Tfv - t)))) + h + 1/(2*alpha)*log(1 + alpha/g*v0^2)
    else
        return -1/alpha*log(abs(cosh(sqrt(alpha*g)*(Tfv - t)))) + h + 1/(2*alpha)*log(1 + alpha/g*v0^2)
    end
end

# the ball velocity when it hits the floor
function vfall_vf(v0, g, h, alpha)
    h_m = x_vf(v0, g, h, alpha, T_vf(v0, g, alpha))
    return sqrt(alpha/g)*sqrt(1 - exp(-2alpha*h_m))
end

### display a comparison 
function simplediff(v0, g, h, alpha, beta, N, T, kind)
    
    # time initialization
    if(T < 0)
        T = T_vf(v0, g, alpha)
    end
    t = [i*T/(N-1) for i=0:(N-1)]

    # selection
    if(kind == "speed")
        # construction of the velocities
        Vsf = [v_sf(v0, g, t[i]) for i=1:N]
        Vff = [v_ff(v0, g, beta, t[i]) for i=1:N]
        Vvf = [v_vf(v0, g, alpha, t[i]) for i=1:N]

        # displaying canvas
        plot(t, [Vsf, Vff, Vvf], grid=false, w=[1 1 1], color = ["blue" "green" "red"], linestyle=[:solid :solid :solid], label=["vsf" "vff" "vvf"])
    elseif(kind == "position")
        # construction of the velocities
        Xsf = [x_sf(v0, g, h, t[i]) for i=1:N]
        Xff = [x_ff(v0, g, h, beta, t[i]) for i=1:N]
        Xvf = [x_vf(v0, g, h, alpha, t[i]) for i=1:N]

        # displaying canvas
        plot(t, [Xsf, Xff, Xvf], ylim=[0.0, max(maxi(Xsf), maxi(Xff), maxi(Xvf))], grid=false, w=[1 1 1], color = ["blue" "green" "red"], linestyle=[:solid :solid :solid], label=["xsf" "xff" "xvf"])
    end
end

# let's see what we get in the end with the 3 different models with an example
begin
    v0 = 5.0
    g = 9.81
    h = 3.3
    alpha = 10.0
    beta = 9.5
    N = 100
    T = 1.2
    simplediff(v0, g, h, alpha, beta, N, T, "position")
end