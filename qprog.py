from quantum import *

# print(CNOTR(controlbitindex=1))
# print(CNOTR(controlbitindex=1).span)

# q = qprogram(2, name = 'refactor testing')

# q.addgates(0, [HGATE(), HGATE()])
# q.addgates(1, [IGATE(), CNOTR()])

# q.compile()

# q.run(graph = True, terminal = True)

bruh = Matrix(3, 3)
bruh.rows[0][0] = comp(1)
bruh.rows[1][0] = comp(1)
bruh.rows[2][0] = comp(1)

bruh.rows[0][1] = comp(2)
bruh.rows[1][1] = comp(3)
bruh.rows[2][1] = comp(4)

bruh.rows[0][2] = comp(4)
bruh.rows[1][2] = comp(9)
bruh.rows[2][2] = comp(16)

# a = Matrix(2, 2)
# a.rows[0][0] = comp(1)
# a.rows[0][1] = comp(2)
# a.rows[1][0] = comp(3)
# a.rows[1][1] = comp(4)
print(bruh)
print()
print(bruh.inverse())