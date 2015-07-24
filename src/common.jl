# common facilities

### value type conversion

f64(x::Real) = Float64(x)


### interface helper

function _req_method(fex)
    fex.head == :call ||
        error("Macro @req_method should be applied to a call expression.")

    f = fex.args[1]

    n = length(fex.args) - 1
    argsyms = Array(Symbol, n)
    for i = 1:n
        aex = fex.args[i+1]
        if isa(aex, Symbol)
            argsyms[i] = aex
        elseif aex.head == :(::)
            argsyms[i] = aex.args[1]
        else
            error("Invalid function argument $(aex).")
        end
    end
    argtup = Expr(:tuple, argsyms...)
    :($fex = throw(MethodError($f, $argtup)))
end

macro req_method(fex)
    esc(_req_method(fex))
end
