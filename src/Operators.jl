# operator library

module Operators

using LinearAlgebra

export Operator, Gate, Measurement, GATES, PARAM_GATES

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
                             0 0 0 1]
                   )

const PARAM_GATES = Dict(:RX => θ -> exp(im * θ * GATES[:X]),

                         :RY => θ -> exp(im * θ * GATES[:Y]),

                         :RZ => θ -> exp(im * θ * GATES[:Z]),

                         :PHASE => α -> [1 0; 0 exp(im * α)]
                         )

const C = Complex{Float32}

const BASIS = (C.([1, 0]), C.([0, 1]))

abstract type Operator end

mutable struct Measurement <: Operator
    k::Int
    basis::NTuple{2, Vector{C}}

    Measurement(k::Int, basis=BASIS) = new(k, basis)
end

function (M::Measurement)(ψ::Vector{C}, N::Int)
    β₁ = M.basis[1]
    M₁ = β₁ * β₁'
    M̃₁ = tensor(M₁, M.k, N)
    P₁ = real(ψ' * M̃₁ * ψ)
    u = rand(Float32)
    if u < P₁
        p = sqrt(P₁)
        return (M̃₁ / p * ψ, 0)
    else
        β₂ = M.basis[2]
        M₂ = β₂ * β₂'
        M̃₂ = tensor(M₂, M.k, N)
        p = sqrt(1 - P₁)
        return (M̃₂ / p * ψ, 1)
    end
end

struct Gate <: Operator
    k::Tuple{Vararg{Int}}
    rep::Matrix{C}
    name::Symbol
    param::Union{C, Nothing}

    Gate(tag::Tuple{Symbol,Int}) = begin
        name, q = tag;
        new((q,), GATES[name], name, nothing)
    end

    Gate(tag::Tuple{Symbol,Int,Int}) = begin
        name, q1, q2 = tag;
        new((q1, q2), GATES[name], name, nothing)
    end

    Gate(tag::Tuple{Symbol,C,Int}) = begin
        name, p, q = tag;
        new((q,), PARAM_GATES[name](p), name, p)
    end

    Gate(tag::Tuple{Symbol,C,Int,Int}) = begin
        name, p, q1, q2 = tag;
        new((q1, q2), PARAM_GATES[name](p), name, p)
    end
end

(G::Gate)(ψ::Vector{C}, N::Int) = tensor(G.rep, G.k, N) * ψ

# tensor the matrix U of order n up to order N,
# where U acts on qubits: j, j + 1,... k; j < k
#
# and |ψ₀> = |0>₁ ⊗ |0>₂ ⊗ ... ⊗ |0>ₙ
#
# for example, with 2 qubits in Q₁ ⊗ Q₂, let's entangle them:
#
# let U₁ = CX(2, 1)
#        = σ(2, 1)CX(1, 2)σ(1, 2)
#
#     U₁(|00> + |10>) = 
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

# Uₖ -> I₁ ⊗ ... ⊗ I(k - 1 times) ⊗ ... ⊗ I(k - 1)
#       ⊗ Uₖ ⊗ I(k + n - 1) ⊗ ... ⊗ I(N - k - n + 1 times) ⊗ ... ⊗ I(N)
#
function tensor(U::Matrix{C}, k::Int, N::Int)
    n = Int(log(2, size(U, 1)))
    L = diagm(ones(C, 2^(k - 1)))
    R = diagm(ones(C, 2^(N - k - n + 1)))
    L ⊗ C.(U) ⊗ R
end

# Qᵢ -> Q[i + 1]
τ(i, N) = tensor(C.(GATES[:SWAP]), i, N)

# j > k => Qₖ -> Qⱼ
σ(j, k, N) = k < j ? *([τ(k + j - i - 1, N) for i = k:j-1]...) : C(1)

# j < k => Qⱼ -> Q(k - 1)
σ′(j, k, N) = σ(k - 1, j, N)

⊗(x, y) = kron(x, y)

end
