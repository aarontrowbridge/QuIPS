push!(LOAD_PATH, pwd())

using QuIPS

⊗(x, y) = kron(x, y)

# ψ₀ = H * [0,1] ⊗ [1,0] ⊗ [1,0]

quip = [
    # entangeling alice and bob's qubits
    (:H, 2),
    (:CX, (2, 3)),

    # preparing qubit to be teleported
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

println()

M₁ = Int(QC.out[1])
M₂ = Int(QC.out[2])

X̂(i) = Gate((:X, i))
Ẑ(i) = Gate((:Z, i))

X = GATES[:X]
H = GATES[:H]

ψ′ = (X^M₁ * [1,0]) ⊗ (X^M₂ * [1,0]) ⊗ (H * X * [1,0])

operate!(QC, X̂(3)^M₁)
operate!(QC, Ẑ(3)^M₂)

for d in QC.wfn println(d) end

println()

for d in ψ′ println(d) end



