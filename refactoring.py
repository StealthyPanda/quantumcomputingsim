from quantum import *

# mat = Matrix(3, 3)

# mat.rows[0][0] = comp(1)
# mat.rows[1][0] = comp(1)
# mat.rows[2][0] = comp(1)


# mat.rows[0][1] = comp(2)
# mat.rows[1][1] = comp(3)
# mat.rows[2][1] = comp(4)


# mat.rows[0][2] = comp(4)
# mat.rows[1][2] = comp(9)
# mat.rows[2][2] = comp(16)

# # a = Matrix(2, 3)
# a = Matrix([
#     [comp(1), comp(2), comp(3)],
#     [comp(4), comp(5), 6]
# ])

# # b = Matrix(3, 2)
b = Matrix([
    [comp(1), 2],
    [comp(3), comp(4)],
    [comp(4), comp(5)]
])

# c = [2, 3]
# d = [comp(0, 1), comp(2)]

# m = mtensor(a, b)
# print(m, m.shape)
# print(a)
# print(type(pow(-1, 0.5)))
# print(MEASURE([pow(2, -0.5), 0, pow(2, -0.5), 0]))
# run(
#     500, [pow(2, -0.5), 0, pow(2, -0.5), 0], graph = True
# )
# print(b)
# print()
# print(mtensor(b, Matrix([[1]])))

qprog = qprogram(3, name='refactor testing')

something = Matrix(8, 8, id = 's')

qprog.addgates(0, [something])
qprog.addgates(1, [igate, cnot1])
qprog.addgates(2, [igate, hgate])

qprog.compile()

qprog.run(graph = True)