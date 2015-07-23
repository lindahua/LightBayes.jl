using Distributions
using BayesBase
using Base.Test


function test_efd(ed::EFUnivariateDistribution, d0::UnivariateDistribution)
    D = typeof(d0)
    dc = convert(D, ed)
    @test dc == d0
    @test convert(Distribution, ed) == d0
    @test convert(typeof(ed), d0) == ed

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
    (EFD.Normal(2.0), Normal(2.0, 1.0)),
    (EFD.Normal(3.0, 4.0), Normal(0.75, 0.5))
    ]

    println("    testing ", ed)
    test_efd(ed, d0)
end
