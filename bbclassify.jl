using Pkg
Pkg.add("QuadGK")
using QuadGK

function etl(mu, sigma, reliability, minimum, maximum)
    ((mu - minimum) * (maximum - mu) - (reliability * sigma)) / (sigma * (1 - reliability))
end

function k(mu, sigma, reliability, length) 
    sigma_e = sigma * (1 - reliability)
    num = length * ((length - 1) * (sigma - sigma_e) - length * sigma + mu * (length - mu))
    den = 2 * (mu * (length - mu) - (sigma - sigma_e))
    num / den
end

function beta(a, b)
    quadgk(x -> x^(a - 1) * (1 - x)^(b - 1), 0, 1)[1]
end

function dbeta(x, a, b, l, u)
    if x < l || x > u
        0
    else
        (1 / beta(a, b))  * (((x - l)^(a - 1) * (u - x)^(b - 1)) / (u - l)^(a + b - 1))
    end
end

function dbinom(p, n, N)
    binomial(N, n) * (p^n * (1 - p)^(N - n))
end

function dcbinom(p, n, N, k)
    a = dbinom(p, n, N)
    b = dbinom(p, n, N - 2)
    c = dbinom(p, n - 1, N - 2)
    d = dbinom(p, n - 2, N - 2)
    e = k * p * (1 - p)
    a - e * (b - 2*c + d)
end
