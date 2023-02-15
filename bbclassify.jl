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

function dfac(x, r)
    x_o = copy(x)
    for i in 1:length(x)
        if r <= 1
            x[i] = x[i]^r
        else
            for j in 1:r
                if j > 1
                    x[i] = x[i] * (x_o[i] - j + 1)
                end
            end
        end
    end
    x
end

function tsm(x, n, k)
    m = [0.0, 0.0, 0.0, 0.0]
    for i in 1:4
        if i == 1
            m[i] = mean(x) / n
        else
            y = copy(x)
            a = mean(dfac(y, i))[1]
            b = dfac([n - 2], i - 2)[1]
            c = 1 / dfac([n], 2)[1]
            m[i] = (a / b) * c
        end
    end
    m
end

function bbintegrate1(a, b, l, u, N, n, k, lower, upper; method = "ll")
    if method == "ll"
        quadgk(x -> dbeta(x, a, b, l, u) * dbinom(x, n, N), lower, upper)[1]
    else
        quadgk(x -> dbeta(x, a, b, l, u) * dcbinom(x, n, N, k), lower, upper)[1]
    end
end

function bbintegrate1(a, b, l, u, N, n1, n2, k, lower, upper; method = "ll")
    if method == "ll"
        quadgk(x -> dbeta(x, a, b, l, u) * dbinom(x, n1, N) * dbinom(x, n2, N), lower, upper)[1]
    else
        quadgk(x -> dbeta(x, a, b, l, u) * dcbinom(x, n1, N, k) * dcbinom(x, n2, N, k), lower, upper)[1]
    end
end

function betaparameters(x, n, k, model, l, u)
    m = tsm(x, n, k)
    s2 = m[2] - m[1]^2
    g3 = (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3) / (math.sqrt(s2)^3)
    g4 = (m[4] - 4 * m[1] * m[3] + 6 * m[1]^2 * m[2] - 3 * m[1]^4) / (math.sqrt(s2)**4)
    if model == 4
        r = 6 * (g4 - g3^2 - 1) / (6 + 3 * g3^2 - 2 * g4)
        if g3 < 0
            a = r / 2 * (1 + (1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1))))^0.5)
            b = r / 2 * (1 - (1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1))))^0.5)
        else
            b = r / 2 * (1 + (1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1))))^0.5)
            a = r / 2 * (1 - (1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1))))^0.5)
        l = m[1] - ((a * (s2 * (a + b + 1))^0.5) / (a * b)^0.5)
        u = m[1] + ((b * (s2 * (a + b + 1))^0.5) / (a * b)^0.5)
        end
    end
    if model == 2
        a = ((l - m[1]) * (l * (m[1] - u) - m[1]^2 + m[1] * u - s2)) / (s2 * (l - u))
        b = ((m[1] - u) * (l * (u - m[1]) + m[1]^2 - m[1] * u + s2)) / (s2 * (u - l))
    end
    [a, b, l, u]
end
