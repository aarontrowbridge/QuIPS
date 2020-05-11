# Quantum Virtual Machine Library

module QVM

export Circuit, run!, step!, reset!, operate!

using Operators

const C = Complex{Float32}

mutable struct Circuit
    wfn::Vector{C}
    ops::Vector{Operator}
    out::BitVector
    pos::Int
    N::Int

    function Circuit(quip::Vector, N::Int)
        wfn = zeros(2^N); wfn[1] = 1
        ops = compile(quip)
        out = BitVector(undef, N)
        new(wfn, ops, out, 0, N)
    end
end

function compile(quip::Vector)
    Vs = []
    for tag in quip
        if tag[1] == :MEASURE
            push!(Vs, Measurement(tag[2]))
        else
            push!(Vs, Gate(tag))
        end
    end
    Vs
end

function run!(C::Circuit)
    for V in C.ops
        C.pos += 1
        evolve!(C, V)
    end
end

function step!(C::Circuit)
    if C.pos < length(C.ops)
        C.pos += 1
        V = C.ops[C.pos]
        evolve!(C, V)
    end
end

function reset!(C::Circuit)
    C.wfn = zero.(C.wfn); C.wfn[1] = 1
    C.out = BitVector(undef, C.N)
    C.pos = 0
end

function evolve!(C::Circuit, V::Operator)
    if typeof(V) == Measurement
        ψ, outcome = V(C.wfn, C.N)
        C.out[V.k] = outcome
        C.wfn .= ψ
    else
        ψ = V(C.wfn, C.N)
        C.wfn .= ψ
    end
end

function operate!(C::Circuit, V::Operator)
    insert!(C.ops, C.pos + 1, V)
    step!(C)
end

end
