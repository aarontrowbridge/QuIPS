#
# script, or "quip", to run the CIRQ
#

push!(LOAD_PATH, pwd())

using QuIPS

M₁ = Measurement(1)

# quip = [
#     (:H, 1)
# ]
# cirq = Circuit(superpose, 1)
# run!(cirq)
# println(abs.(cirq.wfn))
# ψ, out = M₁(cirq.wfn, cirq.N)
# println(out)

quip = [
    (:H, 2),
    (:CNOT, 2, 1),
    (:MEASURE, 1),
    (:MEASURE, 2)
]
cirq = Circuit(entangle, 2)
for i = 1:10
    run!(cirq)
    println(cirq.out)
    reset!(cirq)
end
N = 10000
gs = 0
for i = 1:N
    run!(cirq)
    if cirq.out[1] == 1
        global gs += 1
    end
    reset!(cirq)
end
println(gs/N)



