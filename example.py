from quantum import *


# state = tensor(HAD(qbit(0)), qbit(0))
# state = tensor(state, HAD(qbit(1)))

# print(state, "\n\n")
# print(HAD(state))
print(HGATE(2), '\n\n')
h = HGATE(2)

newone = mtensor(HGATE(1), HGATE(1))
bruh = True

for each in range(4):
    for i in range(4):
        bruh = bruh and (h.rows[each][i] == newone.rows[each][i])

m = Matrix(2, 1)
m.rows[0][0] = comp(1, 0)
m.rows[1][0] = comp(0, 0)

print(mtensor( m , HGATE(1)))