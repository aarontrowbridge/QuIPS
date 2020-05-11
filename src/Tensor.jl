#
# Tensor the matrix V of order n up to order N,
# where V acts on qubits: j, j + 1,... k; j < k
#
# and |ψ₀> = |0>₁ ⊗ |0>₂ ⊗ ... ⊗ |0>(N)
#
# V(k) -> I₁ ⊗ ... ⊗ I(k - 1 times) ⊗ ... ⊗ I(k - 1)
#       ⊗ V(k) ⊗ I(k + n - 1) ⊗ ... ⊗ I(N - k - n + 1 times) ⊗ ... ⊗ I(N)
#

module Tensor

using LinearAlgebra

export tensor

const C = Complex{Float32}

function tensor(V::Matrix{C}, indx::Vector{Int}, N::Int)
    if length(indx) == 2
        j, k = indx[1], indx[2]
        if j > k
            # Qₖ ⊗ Qⱼ -> Qⱼ ⊗ Q(j - 1)
            Ṽ = tensor(V, j, N)
            F, B = σ(j, k, N)
            return B * Ṽ * F
        else
            # Qⱼ ⊗ Qₖ -> Q(k - 1) ⊗ Qₖ
            Ṽ = tensor(V, k, N)
            F, B = σ(k - 1, j, N)
            return B * Ṽ * F
        end
    # else
    #     i, j, k = indx[1], indx[2], indx[3]
    #     if i > j && j > k
    #         Ṽ = tensor(V, j - 1, N)
    end
end


function tensor(V::Matrix{C}, k::Int, N::Int)
    n = Int(log(2, size(V, 1)))
    L = diagm(ones(C, 2^(k - n)))
    R = diagm(ones(C, 2^(N - k)))
    L ⊗ V ⊗ R
end

⊗(x, y) = kron(x, y)

const SWAP = [1 0 0 0;
              0 0 1 0;
              0 1 0 0;
              0 0 0 1]

# Q(i) -> Q(i + 1)
τ(i, N) = tensor(C.(SWAP), i + 1, N)

# k < j => Q(k) -> Q(j)
σ(j, k, N) =  begin
    τs = [τ(k + j - i - 1, N) for i = k:j-1];
    reverse_τs = reverse(τs);
    k < j ? (*(τs...), *(reverse_τs...)) : (C(1), C(1))
end

end
