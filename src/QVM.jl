# Quantum Virtual Machine Library

module QVM

export Circuit, run!, step!, reset!

using Operators

mutable struct Circuit
    N::Int
    wfn::Vector{Complex{Float32}}
    ops::Vector{Operator}
    out::BitVector
    pos::Int

    function Circuit(quip, N::Int)
        wfn = zeros(2^N); wfn[1] = 1
        ops = compile(quip)
        out = BitVector(undef, N)
        new(N, wfn, ops, out, 0)
    end
end

function compile(quip)
    oprs = []
    for tag in quip
        if tag[1] == :MEASURE
            push!(oprs, Measurement(tag[2]))
        else
            push!(oprs, Gate(tag))
        end
    end
    oprs
end

function run!(cirq::Circuit)
    for opr in cirq.ops
        if typeof(opr) == Measurement
            ψ, outcome = opr(cirq.wfn, cirq.N)
            cirq.out[opr.k] = outcome
            cirq.wfn .= ψ
        else
            ψ = opr(cirq.wfn, cirq.N)
            cirq.wfn .= ψ
        end
        cirq.pos += 1
    end
end

function step!(cirq::Circuit)
    if cirq.pos < length(cirq.ops)
        cirq.pos += 1
        opr = cirq.ops[cirq.pos]
        if typeof(opr) == Measurement
            ψ, outcome = opr(cirq.wfn, cirq.N)
            cirq.out[opr.k] = outcome
            cirq.wfn .= ψ
        else
            ψ = opr(cirq.wfn, cirq.N)
            cirq.wfn .= ψ
        end
     end
end

function reset!(cirq::Circuit)
    d = length(cirq.wfn)
    cirq.wfn = zeros(Float32, d); cirq.wfn[1] = 1
    cirq.out = BitVector(undef, cirq.N)
    cirq.pos = 0
end

end
