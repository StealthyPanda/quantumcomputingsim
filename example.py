from quantum import *


# create a qbit using the qbit(x) function, where x is the default value of the bit:
q0 = qbit(1)
q1 = qbit(0)
#qubit 1 is represented by [0 1] while qubit 0 is [1 0]


#all qbits are just a list of 2 elements, elements being either real or complex
#complex numbers in this library are handled with `comp` class, and Matrices with `Matrix` class
print("q0:", q0, ", q1:", q1)


#gates are applied on qubits using the corresponding gate funtions.
#as of now, IDEN (identity), NOT (not), HAD (hadamard), CNOT (controlled not) gates have been implemented. 
#however, after application of CNOT gate I still haven't figured out how to use other gates

#to see the outcome of a computation, use the MEASURE function
#Measurement can be done on a single qubit or a set of qubits. If multiple qubits are to be measured,
#then we must pass the inner tensor product of these qubits to MEASURE. tensor product can be 
#calculated using the tensor function.

#example: measures a single qubit (here the qubit is 1, so measurement is always 1, try running this line a few times)
print(MEASURE(q0))
# print(MEASURE(q1))

#example: measuring a set of qubits:
q0andq1 = tensor(q0, q1)
print(MEASURE(q0andq1))

#example: displaying measurements in a readable way
#use the extract function to get the value of a single qubit from a measurement of set of qubits
measurement = MEASURE(q0andq1)
q0m = extract(measurement, 0) #here 0 is the qbitindex as it is the 0th element in the tensor product
q1m = extract(measurement, 1)
print("q0m: ", q0m, ", q0m: ", q1m)