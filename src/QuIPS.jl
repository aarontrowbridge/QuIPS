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
export Control
export Measurement
export GATES
export PARAM_GATES


using Tensor

export tensor
export ⊗

end
