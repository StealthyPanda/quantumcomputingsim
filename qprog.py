from quantum import *

test = qprogram(3)

test.addgates(0, [HAD])
test.addgates(1, [HAD, NOT])
test.addgates(2, [HAD, HAD, NOT])

aq = NOT(HAD(HAD(qbit(0))))

print(extract(MEASURE(test.getstate(0)), 0), extract(MEASURE(test.getstatetensor()), 0))

# for each in range(10):
#     m = test.measure()
#     print(m.index(1), extract(m, 0), extract(m, 1), extract(m, 2), extract(MEASURE(aq), 0), m)

# for each in range(10):
#     print((test.measure(0)), (test.measure(1)), (test.measure(2)))