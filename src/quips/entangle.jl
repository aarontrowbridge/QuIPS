push!(LOAD_PATH, pwd())

using QuIPS

quip = [
    (:H, 1),
    (:CX, [1, 2]),
    (:MEASURE, 1),
    (:MEASURE, 2)
]

QC = QCircuit(quip, 2)

for i = 1:5
    run!(QC)
    println(QC.out)
    reset!(QC)
end

println()

N = 1000
gs = 0
for i = 1:N
    run!(QC)
    if QC.out[1] == 1
        global gs += 1
    end
    reset!(QC)
end

println("% of correlated outcomes = 1 for $N runs: ", gs/N)



