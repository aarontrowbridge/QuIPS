# Quantum Information Processing System Library

push!(LOAD_PATH, pwd())

using Operators

mutable struct QVM
    N::Int
    wfn::Vector{Complex{Float32}}
    ops::Vector{Operator}
    out::BitVector
    pos::Int

    function QVM(quip, N::Int)
        wfn = zeros(2^N); wfn[1] = 1
        ops = compile(quip)
        out = BitVector(undef, N)
        new(N, wfn, ops, out, 0)
    end
end

function compile(quip)
    gates = []
    measurements = []
    for (j, tag) in enumerate(quip)
        if tag[1] == :MEASURE
            push!(measurements, Measurement(tag[2], j))
        else
            push!(gates, Gate(tag, j))
        end
    end
    ops = [gates; measurements]
    sort!(ops, by=(opr -> opr.loc))
    ops
end

function run!(qvm::QVM)
    for opr in qvm.ops
        if typeof(opr) == Measurement
            ψ, outcome = opr(qvm.wfn, qvm.N)
            qvm.out[opr.k] = outcome
            qvm.wfn .= ψ
        else
            ψ = opr(qvm.wfn, qvm.N)
            qvm.wfn .= ψ
        end
        qvm.pos += 1
    end
end

function step!(qvm::QVM)
    if qvm.pos < length(qvm.ops)
        qvm.pos += 1
        opr = qvm.ops[qvm.pos]
        if typeof(opr) == Measurement
            ψ, outcome = opr(qvm.wfn, qvm.N)
            qvm.out[opr.k] = outcome
            qvm.wfn .= ψ
        else
            ψ = opr(qvm.wfn, qvm.N)
            qvm.wfn .= ψ
        end
     end
end

function reset!(qvm::QVM)
    d = length(qvm.wfn)
    qvm.wfn = zeros(Float32, d); qvm.wfn[1] = 1
    qvm.out = BitVector(undef, qvm.N)
    qvm.pos = 0
end

circuit = [(:H, 1),
           (:CNOT, 2),
           (:MEASURE, 1),
           (:MEASURE, 2)]

qvm = QVM(circuit, 2)


M = Measurement(1)

for i = 1:10
    run!(qvm)
    println(qvm.out)
    reset!(qvm)
end


N = 10000
gs = 0
for i = 1:N
    run!(qvm)
    if qvm.out[1] == 1
        global gs += 1
    end
    reset!(qvm)
end
println(gs/N)






