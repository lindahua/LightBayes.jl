tests = ["efd_univariate"]

for t in tests
    fpath = "$t.jl"
    println("running $fpath ...")
    include(fpath)
end
