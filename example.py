from quantum import *



test = qprogram(3)
test.addgates(0, [HGATE(), HGATE(),  IGATE()])
test.addgates(1, [HGATE(), HGATE(), NGATE()])
test.addgates(2, [IGATE(), CNOTGATE(), NGATE()])
# test.addgates(4, [HGATE()])
# test.addgates(5, [HGATE()])


test.compile()
# test.run(verbose = True)
test.run(graph = True)

#victory