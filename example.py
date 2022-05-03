from quantum import *

# print(comp.getcomplex(comp(1, 0)))
# print(FLIP(tensor(qbit(0), qbit(1))))
s = CNOT(SHIFT(HAD(qbit(1)), PI/3), (qbit(0)))

s = SHIFT(s, PI/3)

print(s, MEASURE(s))

run(5, s)