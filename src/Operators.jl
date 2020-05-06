# operator library

module Operators

export Operator, Gate, Measurement

const gates = Dict(:I => [1 0;
                          0 1],

                   :X => [0 1;
                          1 0],

                   :Y => [0 -im;
                          im 0],

                   :Z => [1 0;
                          0 -1],

                   :H => 1/sqrt(2) * [1 1;
                                      1 -1],

                   :CNOT => [1 0 0 0;
                             0 1 0 0;
                             0 0 0 1;
                             0 0 1 0],

                   :SWAP => [1 0 0 0;
                             0 0 1 0;
                             0 1 0 0;
                             0 0 0 1]
                   )

const param_gates = Dict(:RX => θ -> exp(im * θ * gates[:X]),

                         :RY => θ -> exp(im * θ * gates[:Y]),

                         :RZ => θ -> exp(im * θ * gates[:Z]),

                         :PHASE => α -> [1 0; 0 exp(im * α)]
                         )


abstract type Operator end

struct Gate <: Operator
    op::Matrix{Complex{Float32}}
    pos::Int
    name::Symbol
    param::Union{Complex{Float32}, Nothing}
    target::Tuple{Vararg{Int}}

    Gate(opr::Tuple{Symbol,Int}, pos::Int) = begin
        name, q = opr;
        new(gates[name], pos, name, nothing, (q,))
    end

    Gate(opr::Tuple{Symbol,Int,Int}, pos::Int) = begin
        name, q1, q2 = opr;
        new(gates[name], pos, name, nothing, (q1, q2))
    end

    Gate(opr::Tuple{Symbol,Complex{Float32},Int}, pos::Int) = begin
        ((name, p, q), pos) = opr;
        new(param_gates[name](p), pos, name, p, (q,))
    end

    Gate(opr::Tuple{Symbol,Complex{Float32},Int,Int}, pos::Int) = begin
        ((name, p, q1, q2), pos) = opr;
        new(param_gates[name](p), pos, name, p, (q1, q2))
    end
end

(G::Gate)(ket::Vector{Complex{Float32}}) = begin
    n = Int(log(2, length(ket)));
    return tensor(G, n) * ket
end


mutable struct Measurement <: Operator
    pos::Int
    basis::NTuple{2, Vector{Complex{Float32}}}
    target::Int
    outcome::Union{Bool, Nothing}
    orthonormal::Bool

    Measurement(target, pos) = new(pos, ([1,0],[0,1]), target, nothing, true)
end

function(M::Measurement)(ket::Vector{Complex{Float32}})
    n = Int(log(2, length(ket)))
    β1 = M.basis[1]
    β2 = M.basis[2]
    op1 = lift(β1 * β1', (M.target,), n)
    if M.orthonormal
        P1 = sqrt(real(ket' * op1 * ket))
    else
        P1 = sqrt(real(ket' * op1' * op1 * ket))
    end
    u = rand(Float32)
    if u < P1
        M.outcome = false
        return op1 * ket / P1
    else
        op2 = lift(β2 * β2', (M.target,), n)
        P2 = sqrt(1 - P1^2)
        M.outcome = true
        return op2 * ket / P2
    end
end

function tensor(G::Gate, n::Int)
    if length(G.target) == 1 || abs(G.target[2] - G.target[1]) == 1
        return lift(G.op, G.target, n)
    else
        j, k = G.target
        if j > k
            return P_(k, j, n) * lift(G.op, G.target, n) * P_(j, k, n)
        else
            return _P(k, j, n) * lift(G.op, G.target, n) * _P(j, k, n)
        end
    end
end

P_(j, k, n) = *([lift(gates[:SWAP], (j+k-i-2, j+k-i-1), n) for i = k:j-2]...)
_P(j, k, n) = lift(gates[:SWAP], (k-1, k), n) * P(j, k, n)

lift(U::Matrix{Complex{Float32}}, q::Tuple{Vararg{Int}}, n::Int) = begin
    L = nkron(gates[:I], minimum(q) - 1);
    R = nkron(gates[:I], n - maximum(q));
    kron(kron(L, U), R)
end

nkron(M, n) = n < 1 ? 1 : kron(M, nkron(M, n - 1))




end
