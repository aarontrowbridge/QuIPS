# QuIPS
**Quantum Information Processing System**

Bertrand Russel said *"War does not determine who is right, only who is left"*.

We can interpret this as, *Quantum Mechanics does not determine which outcome is right, only which operators are to the left*.

That motto will be the starting point of this project; modeling a quantum computer as a sequence of gate operations on an initital tensor product state of n Qubits.  

Julia is being used because it is better then whatever language you want to compare it to. 

## Structs

### Qubit

### Ket <: AbstractState

### Density <: AbstractState

### Unitary <: Matrix
(p) -> Unitary(p)

### Single Qubit Gate <: AbstractGate
* G(p)::Unitary # unitary 2x2 matrix
* qbit::Int     # qubit index number
* p::Float      # matrix parameter
* d::Bool       # dagger (conjugate transpose the matrix) boolean; defaulted to *false*
* t::Symbol     # tag for name of gate

### Two Qubit Gate <: AbstractGate
* G(p)::Unitary # tensor product of two single gubit gate matrices (4x4 matrix)
* qbit::Int     # qubit index number
* p::Float      # matrix parameter
* d::Bool       # dagger (conjugate transpose the matrix) boolean; defaulted to *false*
* t::Symbol     # tag for name of gate

## Methods
