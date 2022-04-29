from sys import exit
from quantum import *

q0 = HAD(qbit(1))
q1 = qbit(0)
q2 = HAD(qbit(1))


# state = tensor(CNOT(q0, q1), q2)
res = [[0, 0], [0, 0], [0, 0]]
for each in range(1600):
    q0q1cnot = MEASURE(CNOT(q0, q1))
    #q0 = extract(q0q1cnot, 0)
    q1ex = extract(q0q1cnot, 1)
    # print(q0, q1)
    # res[0][0] += q0[0]
    # res[0][1] += q0[1]

    res[1][0] += q1ex[0]
    res[1][1] += q1ex[1]

q0 = MEASURE(HAD(q0))
for each in range(1600):
    res[2][0] += q0[0]
    res[2][1] += q0[1]

print(res)