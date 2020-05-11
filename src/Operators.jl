# operator library
push!(LOAD_PATH, pwd())

module Operators

using LinearAlgebra, Tensor

export Operator, Gate, Control, Measurement


const GATES = Dict(:I => [1 0;
                          0 1],

                   :X => [0 1;
                          1 0],

                   :NOT => [0 1;
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
                          0 exp(im * π / 4)])


const PARAM_GATES = Dict(:RX => γ -> exp(im * γ * GATES[:X]),

                         :RY => β -> exp(im * β * GATES[:Y]),

                         :RZ => α -> exp(im * α * GATES[:Z]),

                         :PHASE => θ -> [1 0; 0 exp(im * θ)])


const BASIS = ([1, 0],
               [0, 1])


const C = Complex{Float32}


#
#
# Operators
#
#

abstract type Operator end

#
# Measurement Operator
#

struct Measurement <: Operator
    name::Symbol
    basis::NTuple{2, Vector{C}}
    k::Int

    Measurement(k::Int, β=BASIS) = new(Symbol("M"), β, k)
end

function (M::Measurement)(ψ::Vector{C}, N::Int)
    k = M.k
    β₀ = M.basis[1]
    M₀ = C.(β₀ * β₀')
    M̃₀ = tensor(M₀, k, N)
    P₀ = abs(ψ' * M̃₀ * ψ)
    u = rand(Float32)
    if u < P₀
        p = sqrt(P₀)
        return M̃₀/p * ψ, 0
    else
        β₁ = M.basis[2]
        M₁ = C.(β₁ * β₁')
        M̃₁ = tensor(M₁, k, N)
        p = sqrt(1 - P₀)
        return M̃₁/p * ψ, 1
    end
end

#
# Single Qubit Gate
#

struct Gate <: Operator
    name::Symbol
    p::Union{Float32, Nothing}
    k::Int
    U::Matrix{C}

    Gate(tag::Tuple{Symbol,Int}) = begin
        name, k = tag;
        new(name, nothing, k, GATES[name])
    end

    Gate(tag::Tuple{Symbol,Float32,Int}) = begin
        name, p, k = tag;
        new(name, p, k, PARAM_GATES[name])
    end
end

(G::Gate)(ψ::Vector{C}, N::Int) = tensor(G.U, G.k, N) * ψ


#
# Control Operator
#

struct Control <: Operator
    name::Symbol
    gate::Gate
    indx::Tuple{Vararg{Int}}

    function Control(tag::Tuple{Symbol,Tuple{Vararg{Int}}})
        name, indx = tag
        U = Symbol(String(filter!(l -> l != 'C', [String(name)...])))
        k = indx[end]
        gate = Gate((U, k))
        new(name, gate, indx)
    end

    function Control(tag::Tuple{Symbol,Float32,Tuple{Vararg{Int}}})
        name, p, indx = tag
        U = Symbol(String(filter!(l -> l != 'C', [String(name)...])))
        k = indx[end]
        gate = Gate((U, p, k))
        new(name, gate, indx)
    end
end

function (CG::Control)(ψ::Vector{C}, N::Int)
    Cnum = length(filter!(l -> l == 'C',[String(CG.name)...]))
    I = GATES[:I]
    U = CG.gate.U
    P₀ = BASIS[1] * BASIS[1]'
    P₁ = BASIS[2] * BASIS[2]'
    V = P₀ ⊗ I + P₁ ⊗ U
    if Cnum > 1
        for i = 2:Cnum
            I′ = diagm(ones(C, 2^i))
            V = P₀ ⊗ I′ + P₁ ⊗ V
        end
    end
    tensor(V, CG.indx, N) * ψ
end

⊗(x, y)= kron(x, y)

end

