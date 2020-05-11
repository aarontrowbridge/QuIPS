# operator library

module Operators

using LinearAlgebra

export Operator, Gate, Measurement, GATES, PGATES


const GATES = Dict(:I => [1 0;
                          0 1],

                   :X => [0 1;
                          1 0],

                   :Y => [0 -im;
                          im 0],

                   :Z => [1 0;
                          0 -1],

                   :H => [1 1;
                          1 -1] / sqrt(2),

                   :S => [1 0;
                          0 im],

                   :T => [1 0;
                          0 exp(im * π / 4)],

                   :CNOT => [1 0 0 0;
                             0 1 0 0;
                             0 0 0 1;
                             0 0 1 0],

                   :SWAP => [1 0 0 0;
                             0 0 1 0;
                             0 1 0 0;
                             0 0 0 1])


const PGATES = Dict(:RX => γ -> exp(im * γ * GATES[:X]),

                    :RY => β -> exp(im * β * GATES[:Y]),

                    :RZ => α -> exp(im * α * GATES[:Z]),

                    :PHASE => θ -> [1 0; 0 exp(im * θ)])


const COMP_BASIS = ([1, 0],
                    [0, 1])


const C = Complex{Float32}

abstract type Operator end

struct Measurement <: Operator
    k::Int
    basis::NTuple{2, Vector{C}}

    Measurement(k::Int, basis=COMP_BASIS) = new(k, basis)
end

function (M::Measurement)(ψ::Vector{C}, N::Int)
    k = M.k
    β₁ = M.basis[1]
    M₁ = C.(β₁ * β₁')
    M̃₁ = tensor(M₁, k, N)
    P₁ = abs(ψ' * M̃₁ * ψ)
    u = rand(Float32)
    if u < P₁
        p = sqrt(P₁)
        return M̃₁/p * ψ, 0
    else
        β₂ = M.basis[2]
        M₂ = C.(β₂ * β₂')
        M̃₂ = tensor(M₂, M.k, N)
        p = sqrt(1 - P₁)
        return M̃₂/p * ψ, 1
    end
end

struct Gate <: Operator
    name::Symbol
    p::Union{Float32, Nothing}
    k::Tuple{Vararg{Int}}
    U::Matrix{C}

    Gate(gate::Tuple{Symbol,Int}) = new(gate[1],
                                        nothing,
                                        (gate[2],),
                                        GATES[gate[1]])

    Gate(gate::Tuple{Symbol,Int,Int}) = new(gate[1],
                                            nothing,
                                            (gate[2], gate[3]),
                                            GATES[gate[1]])

    Gate(gate::Tuple{Symbol,Float32,Int}) = new(gate[1],
                                                gate[2],
                                                (gate[3],),
                                                PGATES[gate[1]](gate[2]))

    Gate(gate::Tuple{Symbol,Float32,Int,Int}) = new(gate[1],
                                                    gate[2],
                                                    (gate[3], gate[4]),
                                                    PGATES[gate[1]](gate[2]))
end

(G::Gate)(ψ::Vector{C}, N::Int) = tensor(G.U, G.k, N) * ψ

struct Control <: Operator
    name::Symbol
    ktup::Tuple{Vararg{Int}}
    gate::Gate
end






# tensor the matrix U of order n up to order N,
# where U acts on qubits: j, j + 1,... k; j < k
#
# and |ψ₀> = |0>₁ ⊗ |0>₂ ⊗ ... ⊗ |0>ₙ
#
function tensor(U::Matrix{C}, ktup::Tuple{Vararg{Int}}, N::Int)
    if length(ktup) == 1
        k, = ktup
        Ũ = tensor(U, k, N)
        return Ũ
    else
        j, k = ktup
        if j > k
            # Qₖ ⊗ Qⱼ -> Qⱼ ⊗ Q(j + 1)
            Ũ = tensor(U, j - 1, N)
            return σ(k, j, N) * Ũ * σ(j, k, N)
        else
            # Qⱼ ⊗ Qₖ -> Q(k - 1) ⊗ Qₖ
            Ũ = tensor(U, k - 1, N)
            return σ′(k, j, N) * Ũ * σ′(j, k, N)
        end
    end
end

# U(k) -> I₁ ⊗ ... ⊗ I(k - 1 times) ⊗ ... ⊗ I(k - 1)
#       ⊗ U(k) ⊗ I(k + n - 1) ⊗ ... ⊗ I(N - k - n + 1 times) ⊗ ... ⊗ I(N)
#
function tensor(U::Matrix{C}, k::Int, N::Int)
    n = Int(log(2, size(U, 1)))
    L = diagm(ones(C, 2^(k - 1)))
    R = diagm(ones(C, 2^(N - k - n + 1)))
    L ⊗ U ⊗ R
end

# Q(i) -> Q(i + 1)
τ(i, N) = tensor(C.(GATES[:SWAP]), i, N)

# j > k => Q(k) -> Q(j)
σ(j, k, N) = k < j ? *([τ(k + j - i - 1, N) for i = k:j-1]...) : C(1)

# j < k => Q(j) -> Q(k - 1)
σ′(j, k, N) = σ(k - 1, j, N)

⊗(x, y) = kron(x, y)

end
