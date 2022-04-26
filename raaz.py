import sys
from quantum import *



# print(CNOT(qbit(1), qbit(1)))

# sys.exit(0)

# def mod(x) : return x if x >= 0 else (-1 * x)

ten = CNOT(HAD(qbit(1)), HAD(qbit(1)))
# q0 = qbit(0)
# q1 = qbit(0)

# q0 = HAD(q0)

# state = tensor(CNOT(q0, q1), q0)




# ten = state
# print(ten)
# ten = (qbit(0))

shots = 1600 * 4

res = [0 for i in range(len(ten))]

for each in range(shots):
    measurement = MEASURE(ten)
    # measurement, qc = MEASURE(ten), MEASURE(q0)
    res[measurement.index(1)] += 1
    # print("q0:", bin(qc.index(1)), "Measurement:", bin(measurement.index(1)))
ret = ""

for each in range(len(res)):
    ret += (f"|Ψ{bin(each)[2:]}> : {float(res[each]/shots) * 100}%, ")

ret = ret[:-2]
print(ret)