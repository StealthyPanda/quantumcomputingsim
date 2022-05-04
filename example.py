from quantum import *



test = qprogram(3)
test.addgates(0, [HGATE(), HGATE(),  IGATE()])
test.addgates(1, [HGATE(), HGATE(), NGATE()])
test.addgates(2, [IGATE(), CNOTGATE(), NGATE()])


test.compile()
# test.run(verbose = True)
test.run()


#victory