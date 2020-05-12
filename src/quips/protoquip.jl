#
# script, or "quip", to run the CIRQ
#

push!(LOAD_PATH, homedir()*"/Projects/QuIPS/src")

using QuIPS

M₁ = Measurement(1)

quip = [
    (:X, 1),
    (:H, 1)
]

QC = QCircuit(quip, 1)

run!(QC)

println(abs.(QC.wfn))

ψ, out = M₁(QC.wfn, QC.N)

println("Measurement Result: ", out)

