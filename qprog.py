from quantum import *

test = qprogram(3)

test.addgates(0, [IDEN, HAD, IDEN])
test.addgates(1, [IDEN, HAD, NOT])
test.addgates(2, [HAD, HAD, NOT])

aq = NOT(HAD(HAD(qbit(0))))

for each in range(10):
    m = test.measure()
    print(m.index(1), extract(m, 0), extract(m, 1), extract(m, 2), extract(MEASURE(aq), 0), m)