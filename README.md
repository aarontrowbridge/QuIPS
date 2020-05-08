# QuIPS
**Quantum Information Processing System**

Bertrand Russel said

>*"War does not determine who is right, only who is left"*.

We can interpret this as, 

>*Quantum Mechanics does not determine which outcome is right, only which operators are to the left*.

That motto will be the starting point of this project; modeling a quantum computer as a sequence of gate operations on an initital tensor product state of n Qubits.

Julia is being used because it is better then whatever language you want to compare it to. 

## Usage

A quantum compuatation begins with specifying a list of successive gates and measurements.  We shall call this quantum information process a *quip*. Additionally, the number of qubits, *N* must be specified. 

As an example of what a *quip* currently looks like, here is a Julia vector containing the instructions to put 2 qubits into an entangled Bell singlet state and measure them, the results should be correlated.

```
quip = [
    (:H, 2),
    (:CNOT, 2, 1),
    (:MEASURE, 1),
    (:MEASURE, 2)
]
```

A *quip* is currently formatted as a Julia vector of tuples, where each tuple contains a symbol specifying the desired operation, optionally a parameter for parameterized gates and finally the index(es) of the targeted qubit(s). 

The goal for the future is to use macros to read in a quip in a less noisy format.
## To-Do

* controlled operators
* macros
* visualisations
* algorithms 

