# Quantum Information Processing System Library

push!(LOAD_PATH, pwd())

using Operators

const QuIP = Vector{Union{Tuple{Symbol,Vararg{Int}},
                          Tuple{Symbol,Float32,Vararg{Int}}}}

mutable struct QVM
    wfn::Vector{Complex{Float32}}
    ops::Vector{Operator}
    out::BitVector
    pos::Int

    function QVM(quip, n::Int)
        wfn = zeros(Complex{Float32}, 2^n); wfn[1] = 1
        ops = compile(quip)
        out = BitVector(undef, n)
        new(wfn, ops, out, 0)
    end
end

function compile(quip)
    gates = []
    measurements = []
    for (j, opr) in enumerate(quip)
        if opr[1] == :MEASURE
            push!(measurements, Measurement(opr[2], j))
        else
            push!(gates, Gate(opr, j))
        end
    end
    ops = [gates; measurements]
    sort!(ops, by=(op -> op.pos))
    ops
end

function run!(qvm::QVM)
    for op in qvm.ops
        qvm.wfn .= op(qvm.wfn)
        if typeof(op) == Measurement
            qvm.out[op.target] = op.outcome
        end
        qvm.pos += 1
    end
end

function reset!(qvm::QVM)
    n = length(qvm.wfn)
    qvm.wfn = zeros(Float32, n)
    qvm.wfn[1] = 1
    qvm.out = BitVector(undef, Int(log(2, n)))
    qvm.pos = 0
end


circuit = [(:H, 1),
           (:CNOT, 1, 2),
           (:MEASURE, 1)]

QC = QVM(circuit, 2)


run!(QC)
println(QC.wfn)




