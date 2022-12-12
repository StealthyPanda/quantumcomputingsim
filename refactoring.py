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

a = Matrix(2, 3)
a.rows = [
    [comp(1), comp(2), comp(3)],
    [comp(4), comp(5), comp(6)]
]

# b = Matrix(3, 2)
# b.rows = [
#     [comp(1), comp(2)],
#     [comp(3), comp(4)],
#     [comp(4), comp(5)]
# ]

c = [2, 3, 4]


print(a * c)