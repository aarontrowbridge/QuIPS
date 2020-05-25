push!(LOAD_PATH, homedir()*"/Projects/QuIPS/src")

using QuIPS

# ψ₀ = Ĥ₁ * [0,1] ⊗ [1,0] ⊗ [1,0]

quip = [
    # entangeling alice and bob's qubits
    (:H, 2),
    (:CX, (2, 3)),

    # preparing qubit to be teleported by alice to bob
    (:X, 1),
    (:H, 1),

    # teleportation algorithm
    (:CX, (1, 2)),
    (:H, 1),
    (:M, 1),
    (:M, 2)
]

QC = QCircuit(quip, 3)

QC |> run!

m₁ = Int(QC.out[1])
m₂ = Int(QC.out[2])

X̂(i) = Gate((:X, i))
Ẑ(i) = Gate((:Z, i))

operate!(QC, X̂(3)^m₂)
X̂′ = QC.ops[end]
print("  "); show(X̂′)

operate!(QC, Ẑ(3)^m₁)
Ẑ′ = QC.ops[end]
print("  "); show(Ẑ′)

println("\nQC.wfn:\n")
for i in QC.wfn println("  ", i) end

X = Gates().Us[:X]
H = Gates().Us[:H]

ψ′ = (X^m₁ * [1,0]) ⊗ (X^m₂ * [1,0]) ⊗ (H * X * [1,0])

const C = Complex{Float32}

println("\nψ′:\n")
for i in C.(ψ′) println("  ", i) end

M̂₃ = Measurement(3)

operate!(QC, M̂₃)

println("\nψ′′:\n")
for i in QC.wfn println("  ", i) end

println("\nm₃ = ", QC.out[3], "\n")



