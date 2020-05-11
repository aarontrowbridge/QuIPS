push!(LOAD_PATH, pwd())

using QuIPS

function show(V::Operator)
    if typeof(V) == Control
        println("  ", V.name, " ", ["$j " for j in V.indx]...)
    elseif typeof(V) == Gate
        if V.name == :I
            println("  ", V.name)
        else
            println("  ", V.name, " ", V.k)
        end
    else
        println("  MEASURE ", V.k)
    end
end

# ψ₀ = H * [0,1] ⊗ [1,0] ⊗ [1,0]

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

M₁ = Int(QC.out[1])
M₂ = Int(QC.out[2])

X̂(i) = Gate((:X, i))
Ẑ(i) = Gate((:Z, i))

operate!(QC, X̂(3)^M₂)
X̂′ = QC.ops[end]
show(X̂′)

operate!(QC, Ẑ(3)^M₁)
Ẑ′ = QC.ops[end]
show(Ẑ′)

println("\nQC.wfn:")
for i in QC.wfn println("  ", i) end

println()

X = GATES[:X]
H = GATES[:H]

ψ′ = (X^M₁ * [1,0]) ⊗ (X^M₂ * [1,0]) ⊗ (H * X * [1,0])

const C = Complex{Float32}

println("ψ′:")
for i in C.(ψ′) println("  ", i) end

println()
println(" QC.wfn == ψ′ : ", isapprox(QC.wfn, C.(ψ′)))
println()



