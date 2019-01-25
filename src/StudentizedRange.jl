# Provided by "Doomphoenix Qxz" on discourse.julialang.com

using SpecialFunctions: gamma, erf
using ForwardDiff: derivative
#using Plots # needed for running tests, uncomment if you want that
using QuadGK: quadgk

struct StudentizedRange{T<:Real}
    k::T
    ν::T
end

function StudentizedRange{T}(k, ν) where T
    (k, ν) = (Float64(k), Float64(ν))
    return StudentizedRange(k, ν)
end

StudentizedRange(k::T, ν::T) where {T<:Real} = StudentizedRange{T}(k, ν)
StudentizedRange(k::Real, ν::Real) = StudentizedRange(promote(k, ν)...)
StudentizedRange(k::Integer, ν::Integer) = StudentizedRange(Float64(k), Float64(ν))
StudentizedRange(a) = StudentizedRange(a,a)
StudentizedRange() = StudentizedRange(2,2)


function 𝚽(x)
    return (1+erf(x / √2)) / 2
end

function ϕ(x)
    return derivative(𝚽, x)
end

import Distributions.cdf
function cdf(d::StudentizedRange, q)
    function outer(x)

        function inner(u)
            return ϕ(u) * (𝚽(u) - 𝚽(u - q*x))^(d.k-1)
        end
        inner_part = quadgk(inner, -Inf, Inf)[1]
        return inner_part * x^(d.ν-1) * exp(-x^2*d.ν / 2)
    end
    integral = quadgk(outer, 0.0, Inf)[1]
    return integral * (d.k * d.ν^(d.ν/2)) / (gamma(d.ν/2) * 2^(d.ν/2 - 1))

end

import Distributions.pdf
function pdf(d::StudentizedRange, q)

    function outer(x)

        function inner(u)
            return ϕ(u) * ϕ(u - q*x) * (𝚽(u) - 𝚽(u - q*x))^(d.k-2)
        end
        inner_part = quadgk(inner, -Inf, Inf)[1]
        return inner_part * x^d.ν * ϕ(x*√(d.ν))
    end
    integral = quadgk(outer, 0.0, Inf)[1]
    return integral * (√(2π) * d.k * (d.k-1) * d.ν^(d.ν/2)) / (gamma(d.ν/2) * 2^(d.ν/2 - 1))
end

import Distributions.logpdf
logpdf(d::StudentizedRange, q) = log(pdf(d, q))

# To get quantile to work correctly I had to implement my quick naive version of
# the bisection method. I'm not sure why, but trying to do this with Roots.jl
# was WAY too slow. Like 30+ seconds.
function simple_bisection(f::Function, brackets, abstol=10.0^-6, maxeval=1e3)
    if brackets[1] > brackets[2]
        xmax, xmin = brackets
    else
        xmin, xmax = brackets
    end
    @assert f(xmin) * f(xmax) < 0

    a = xmin
    b = xmax
    error = 1
    numeval = 0
    while error > abstol

        numeval += 1
        if numeval > maxeval break end

        fnew = f((a+b)/2)
        error = abs(fnew)
        if fnew * f(a) < 0
            b = (a+b)/2
        elseif fnew * f(b) < 0
            a = (a+b)/2
        else
            throw(BoundsError("Algorithm failed to converge. Sorry :("))
        end

    end
    return (a+b)/2
end


function quantile(d::StudentizedRange, x, tol=1E-6)
    @assert 0.0 <= x < 1.0
    if x == 0.0 return 0.0 end

    return simple_bisection(y -> cdf(d, y) - x, [0.0, 100.0], tol)
end


function runtests()

    plot(legend=:topleft)

    function test_cdf(k, v)
        rangedist = StudentizedRange(k, v)
        xsp = 0.0:0.01:5.0
        cdflist = [cdf(rangedist, xsp[i])[1] for i in 1:length(xsp)]
        plot!(xsp, cdflist, label="k=$k, df=$v")
    end
    test_cdf(2,2)
    test_cdf(2,4)
    test_cdf(2,8)
    test_cdf(3,10)
    test_cdf(10,10)
    test_cdf(10,100)

    title!("Test of Studentized Range Distribution")
    xlabel!("x")
    ylabel!("CDF (x)")
    savefig("Studentized Range Test2.png")

    plot()
    function test_pdf(k, v)
        rangedist = StudentizedRange(k, v)
        xsp = 0.0:0.01:5.0
        pdflist = [pdf(rangedist, xsp[i])[1] for i in 1:length(xsp)]
        plot!(xsp, pdflist, label="k=$k, df=$v")
    end
    test_pdf(2,2)
    test_pdf(2,4)
    test_pdf(2,8)
    test_pdf(3,10)
    test_pdf(10,10)
    test_pdf(10, 100)
    title!("Test of Studentized Range Distribution")
    xlabel!("x")
    ylabel!("PDF (x)")
    savefig("Studentized Range Test.png")
end
