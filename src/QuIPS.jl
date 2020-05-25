push!(LOAD_PATH, homedir()*"/Projects/QuIPS/src")

#
# QuIPS: Quantum Information Processing System
#

module QuIPS


using QVM

export QCircuit
export run!
export step!
export reset!
export operate!


using Operators

export Operator
export Gate
export ControlGate
export Measurement
export Gates


using Tensor

export tensor
export ⊗

end
