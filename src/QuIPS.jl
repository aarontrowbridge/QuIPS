#
# QuIPS: Quantum Information Processing System
#

module QuIPS

include("QVM.jl")
using .QVM

export QCircuit
export run!
export step!
export reset!
export operate!

include("Operators.jl")
using .Operators

export Operator
export Gate
export ControlGate
export Measurement
export Gates

include("Tensor.jl")
using .Tensor

export tensor
export ⊗

end
