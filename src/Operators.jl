# operator library

module Operators

using LinearAlgebra

export Operator, Gate, Measurement, nkron, GATES

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
    loc::Int
    basis::NTuple{2, Vector{C}}

    Measurement(k::Int, loc::Int, basis=BASIS) = new(k, loc, basis)
    Measurement(k::Int, basis=BASIS) = new(k, 0, basis)
end

function (M::Measurement)(ψ::Vector{C}, N::Int)
    β₁ = M.basis[1]
    M₁ = β₁ * β₁'
    M̃₁ = lift(M₁, M.k, N)
    P₁ = real(ψ' * M̃₁ * ψ)
    u = rand(Float32)
    if u < P₁
        p = sqrt(P₁)
        return (M̃₁ / p * ψ, 0)
    else
        β₂ = M.basis[2]
        M₂ = β₂ * β₂'
        M̃₂ = lift(M₂, M.k, N)
        p = sqrt(1 - P₁)
        return (M̃₂ / p * ψ, 1)
    end
end

struct Gate <: Operator
    k::Tuple{Vararg{Int}}
    rep::Matrix{C}
    loc::Int
    name::Symbol
    param::Union{C, Nothing}

    Gate(tag::Tuple{Symbol,Int}, loc::Int) = begin
        name, q = tag;
        new((q,), GATES[name], loc, name, nothing)
    end

    Gate(tag::Tuple{Symbol,Int,Int}, loc::Int) = begin
        name, q1, q2 = tag;
        new((q1, q2), GATES[name], loc, name, nothing)
    end

    Gate(tag::Tuple{Symbol,C,Int}, loc::Int) = begin
        name, p, q = tag;
        new((q,), PARAM_GATES[name](p), loc, name, p)
    end

    Gate(tag::Tuple{Symbol,C,Int,Int}, loc::Int) = begin
        name, p, q1, q2 = tag;
        new((q1, q2), PARAM_GATES[name](p), loc, name, p)
    end
end

(G::Gate)(ψ::Vector{C}, N::Int) = tensor(G.rep, G.k, N) * ψ

function tensor(U::Matrix{C}, ktup::Tuple{Vararg{Int}}, N::Int)
    if length(ktup) == 1
        k, = ktup
        Ũ = lift(U, k, N)
        return Ũ
    else
        j, k = ktup
        if j > k
            Ũ = lift(U, j, N)
            return σ(k, j, N) * Ũ * σ(j, k, N)
        else
            Ũ = lift(U, k, N)
            return σ′(k, j, N) * Ũ * σ'(j, k, N)
        end
    end
end

function lift(U::Matrix, k::Int, N::Int)
    n = Int(log(2, size(U, 1)))
    L = diagm(ones(C, 2^(k - n)))
    R = diagm(ones(C, 2^(N - k)))
    kron(kron(L, C.(U)), R)
end

τ(i, N) = lift(GATES[:SWAP], i, N)

σ(j, k, N) = k < j-1 ? *([τ(j+k-i-2, N) for i = k:j-2]...) : C(1)

σ′(j, k, N) = τ(k-1, N) * π(k, j, N)

nkron(rep, n) = n < 1 ? 1 : kron(rep, nkron(rep, n-1))

end
