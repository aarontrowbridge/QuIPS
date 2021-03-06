#
# Tensor the matrix V of order n up to order N,
# where V acts on qubits: j, j + 1,... k; j < k
#
# and |ψ₀> = |0>₁ ⊗ |0>₂ ⊗ ... ⊗ |0>(N)
#
# V(k) -> I₁ ⊗ ... ⊗ I(k - 1) ⊗ V(k) ⊗ I(k + n - 1) ⊗ ... ⊗ I(N)
#

module Tensor

using LinearAlgebra

export tensor, ⊗

const C32 = Complex{Float32}

function tensor(V::Matrix{C32}, indx::Tuple{Vararg{Int}}, N::Int)
    if length(indx) == 2
        j, k = indx
        if j > k
            # Qₖ ⊗ Qⱼ -> Qⱼ ⊗ Q(j - 1)
            Ṽ = tensor(V, j, N)
            R, L = σ(j, k, N)
        else
            # Qⱼ ⊗ Qₖ -> Q(k - 1) ⊗ Qₖ
            Ṽ = tensor(V, k, N)
            R, L = σ(k - 1, j, N)
        end
        return L * Ṽ * R
    else
        i, j, k = indx
        if i > j && j > k
            Ṽ = tensor(V, i, N)
            R₁, L₁ = σ(i, k, N)
            R₂, L₂ = σ(i - 1, j - 1, N)
        elseif i < j && j > k
            if i < k
                Ṽ = tensor(V, j, N)
                R₁, L₁ = σ(j, k, N)
                R₂, L₂ = σ(j - 2, i, N)
            else
                Ṽ = tensor(V, j, N)
                R₁, L₁ = σ(j, k, N)
                R₂, L₂ = σ(j - 2, i - 1, N)
            end
        elseif i > j && j < k
            if i < k
                Ṽ = tensor(V, k, N)
                R₁, L₁ = σ(k - 1, j, N)
                R₂, L₂ = σ(j - 2, i - 1, N)
            else
                Ṽ = tensor(V, i, N)
                R₁, L₁ = σ(i, k, N)
                R₂, L₂ = σ(i - 1, j, N)
            end
        else
            Ṽ = tensor(V, k, N)
            R₁, L₁ = σ(k - 1, j, N)
            R₂, L₂ = σ(k - 2, i, N)
        end
        return L₁ * L₂ * Ṽ * R₂ * R₁
    end
end

function tensor(V::Matrix{C32}, k::Int, N::Int)
    n = Int(log(2, size(V, 1)))
    L = diagm(ones(C32, 2^(k - n)))
    R = diagm(ones(C32, 2^(N - k)))
    L ⊗ V ⊗ R
end

⊗(x, y) = kron(x, y)

const SWAP = [1 0 0 0;
              0 0 1 0;
              0 1 0 0;
              0 0 0 1]

# Q(i) -> Q(i + 1)
τ(i, N) = tensor(C32.(SWAP), i + 1, N)

# k < j => Q(k) -> Q(j)
σ(j, k, N) = begin
    if k < j
        τs = [τ(k + j - i - 1, N) for i = k:j-1];
        return (*(τs...), *(reverse(τs)...))
    else
        return (C32(1), C32(1))
    end
end

end
