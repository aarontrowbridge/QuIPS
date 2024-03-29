# QuIPS: Quantum Information Processing System

Bertrand Russel said

>*"War does not determine who is right, only who is left"*.

We can interpret this as, 

>*Quantum Mechanics does not determine which state is right, only which operators are to the left*.

That motto will be the starting point of this project; modeling a quantum computer as a sequence of gate operations on an initital tensor product state of *N* Qubits.

## Installation

To install:

* clone this repository
* in julia REPL of cloned directory:
  * enter pkg mode with `]`
  * run `pkg> activate .` to activate directory environment
  * run `pkg> instantiate` to set up package
  
Now package can be used with `using QuIPS` in julia scripts (quips) or right in REPL.

## Usage

A quantum compuatation begins with specifying a list of successive gates and measurements.  We shall call this quantum information program, or process, or whatever, a *quip*. Additionally, the number of qubits, *N* must be specified. 

As an example of what a *quip* currently looks like, here is a Julia vector containing the instructions to put 2 qubits into an entangled Bell singlet state and measure them, the results should be correlated.

```julia
quip = [
    (:H, 2),
    (:CX, (2, 1)),
    (:MEASURE, 1),
    (:MEASURE, 2)
]
```

Here, `:CX = :CNOT`, both are allowed, but I prefer `:CX` as it flows with the naming convention and differentiates quantum compuation, the X gate, from the classical compuation, the NOT gate. "NOT" is ill defined.

A *quip* is currently formatted as a Julia vector of tuples, where each tuple contains a symbol specifying the desired operation, optionally a parameter for parameterized gates and finally the index(es) of the targeted qubit(s), just the index or an array of of indices for multi qubit gates.

```julia
random_quip = [
    (:RX, 3.14, 1),
    (:Y, 2),
    (:CCZ, (3, 2, 4)),
    (:H, 1),
    (:CY, (3, 2))
    (:CCH, (2, 1, 3))
    (:PHASE, 0.5*3.14, 3)
]
```

Once we have our *quip* specified, we need to compile it, and the number of qubits we want to work with, into a Quantum Circuit, for our QVM to work with.  Here we will just run the quantum computer and that looks like this.

```julia
QC = QCircuit(quip, 2)
run!(QC)
```

We can then access information stored in `QC` like `QC.wfn`.

## To-Do

### general
* jupyter lab 
* google colab

### syntax
* macros
* visualisations
* run! functionality

### algorithms 
* ising
* QFFT
* schur transform

### formalisms
* density operator
  * artificial noise
* path integral
  
### optimizations
* GPU
* graph state
  * Gottesman-Knill theorem
* clifford group


