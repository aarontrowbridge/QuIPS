push!(LOAD_PATH, homedir()*"/Projects/QuIPS/src")

#
# QuIPS: Quantum Information Processing System
#

module QuIPS


using QVM

export QCircuit, run!, step!, reset!, operate!


using Operators

export Operator, Gate, Control, Measurement
export GATES, PARAM_GATES


using Tensor

export tensor, ⊗

end
