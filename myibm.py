from sys import exit
from quantum import *

# def bruh(x):
#     return x

# print(bruh(0, 9))

# exit()
state = tensor(CNOT(HAD(qbit(0)), qbit(1)), HAD(qbit(0)))

# run(1600, state)
for each in range(16):
    m = MEASURE(state)
    m2 = MEASURE(m)
    print(m, m2==m, extract(m, 0))