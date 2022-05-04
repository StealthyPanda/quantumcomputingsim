from quantum import *

test = mtensor(IGATE(1), NGATE(1)) ** tensor(qbit(0), qbit(0))

print (tensor(qbit(0), qbit(0)), test)

# state = CNOT(HAD(qbit(0)), qbit(0))

# run(1600, state)

# # state = HAD(state)

# one = Matrix(2, 2)
# one.rows[0][1] = 0
# one.rows[1][0] = 0

# hadonfirstnothingonsecond = mtensor(one, HGATE(1))

# state = hadonfirstnothingonsecond ** state

# run(1600, state)