using Distributions
using BayesModels
using Base.Test

function test_efd(ed::EFUnivariateDistribution, d0::UnivariateDistribution)
    D = typeof(d0)
    ED = typeof(ed)

    @test convert(D, ed) == d0
    @test convert(ED, d0) == ed
    @test convert(Distribution, ed) == d0
    @test convert(EFDistribution, d0) == ed

    n = 100
    x = rand(d0, n)

    lp0 = logpdf(d0, x)
    lp = zeros(n)
    for i = 1:n
        lp[i] = logpdf(ed, x[i])
    end
    @test_approx_eq lp lp0

    @test_approx_eq logpdf(ed, x) lp0
    @test_approx_eq logupdf(ed, x) - logpartition(ed) lp0
end


for (ed, d0) in [
    # Bernoulli
    (EFBernoulli(), Bernoulli(0.5)),
    (EFBernoulli(1.0), Bernoulli(0.7310585786300049)),
    (EFBernoulli(-1.0), Bernoulli(0.2689414213699951)),
    # Beta
    (EFBeta(), Beta(1.0, 1.0)),
    (EFBeta(1.0, 2.0), Beta(2.0, 3.0)),
    # Normal
    (EFNormal(2.0), Normal(2.0, 1.0)),
    (EFNormal(3.0, 4.0), Normal(0.75, 0.5)),
    ]

    println("    testing ", ed)
    test_efd(ed, d0)
end
