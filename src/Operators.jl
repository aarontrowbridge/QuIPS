push!(LOAD_PATH, homedir()*"/Projects/QuIPS/src")

# operator library

module Operators

using LinearAlgebra, Tensor

export Operator, Gate, ControlGate, Measurement
export Gates, addgate, addgates

#
#
# Constants
#
#

const GATES = Dict(:I => [1 0;
                          0 1],

                   :X => [0 1;
                          1 0], :NOT => [0 1;
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

#
# Gates: Set of all usable gates
#

struct Gates
    Us::Dict{Symbol,Matrix}
    PUs::Dict{Symbol, Function} 

    Gates() = new(GATES, PARAM_GATES)
end

#### Functions to add additional gates to the set

# Add one gate
function addgate!(gates:: Gates, U::Tuple{Symbol, Matrix}):
    merge!(gates.Us, Dict(U))
end

function addgate!(gates:: Gates, PU::Tuple{Symbol, Function}):
    merge!(gates.PUs, Dict(PU))
end

# Add a list of gates
function addgates!(gates::Gates, Us::Vector{Tuple{Symbol, Matrix}}):
    for U in Us:
        adgate!(gates, U)
    end
end

function addgates!(gates::Gates, PUs::Vector{Tuple{Symbol, Function}}):
    for PU in PUs:
        adgate!(gates, PU)
    end
end


const COMP_BASIS = ([1, 0],
                    [0, 1])


const C32 = Complex{Float32}


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
    k::Int
    basis::NTuple{2, Vector{C32}}

    Measurement(k::Int, β=COMP_BASIS) = new(k, β)
end

function (M::Measurement)(ψ::Vector{C32}, N::Int)
    k = M.k
    β₀ = M.basis[1]
    M₀ = C32.(β₀ * β₀')
    M̃₀ = tensor(M₀, k, N)
    P₀ = abs(ψ' * M̃₀ * ψ)
    u = rand(Float32)
    if u < P₀
        p = sqrt(P₀)
        return M̃₀/p * ψ, 0
    else
        β₁ = M.basis[2]
        M₁ = C32.(β₁ * β₁')
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
    p::Union{Float32, Nothing} # parameter
    k::Int                     # target qubit
    U::Matrix{C32}

    Gate(tag::Tuple{Symbol,Int}) = begin
        name, k = tag;
        new(name, nothing, k, GATES[name])
    end

    Gate(tag::Tuple{Symbol,Float32,Int}) = begin
        name, p, k = tag;
        new(name, p, k, PARAM_GATES[name](p))
    end
end

(G::Gate)(ψ::Vector{C32}, N::Int) = tensor(G.U, G.k, N) * ψ


#
# ControlGate Operator
#

struct ControlGate <: Operator
    name::Symbol
    gate::Gate
    indx::Tuple{Vararg{Int}}

    function ControlGate(tag::Tuple{Symbol,Tuple{Vararg{Int}}})
        name, indx = tag
        U = Symbol(String(filter!(l -> l != 'C', [String(name)...])))
        k = indx[end]
        gate = Gate((U, k))
        new(name, gate, indx)
    end

    function ControlGate(tag::Tuple{Symbol,Float32,Tuple{Vararg{Int}}})
        name, p, indx = tag
        U = Symbol(String(filter!(l -> l != 'C', [String(name)...])))
        k = indx[end]
        gate = Gate((U, p, k))
        new(name, gate, indx)
    end
end

function (CG::ControlGate)(ψ::Vector{C32}, N::Int)
    I = GATES[:I]
    P̂₀ = COMP_BASIS[1] * COMP_BASIS[1]'
    P̂₁ = COMP_BASIS[2] * COMP_BASIS[2]'
    Û = CG.gate.U
    CÛ = P̂₀ ⊗ I + P̂₁ ⊗ Û
    Cs = length(CG.indx) - 1
    if Cs > 1
        for i = 2:Cs
            Iᵢ = diagm(ones(C32, 2^i))
            CU = P̂₀ ⊗ Iᵢ + P̂₁ ⊗ CÛ
        end
    end
    tensor(CÛ, CG.indx, N) * ψ
end

Î = Gate((:I, 1))

Base.:^(V̂::Operator, n::Int) = n < 1 ? [Î] : [V̂ for i = 1:n]

function Base.show(V::Operator)
    if typeof(V) == ControlGate
        println(V.name, " ", ["$j " for j in V.indx]...)
    elseif typeof(V) == Gate
        if V.name != :I
            println(V.name, " ", V.k)
        end
    else
        println("MEASURE ", V.k)
    end
end

end

