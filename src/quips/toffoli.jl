push!(LOAD_PATH, pwd())

using QuIPS

quip = [
    (:X, 3),
    (:X, 2),
    (:CCX, (3, 2, 1)),
]

QC = QCircuit(quip, 3)

@time run!(QC)

println(abs2.(QC.wfn))
