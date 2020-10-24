from sympy import *

f1 = Symbol('f1')
f2 = Symbol('f2')
f3 = Symbol('f3')
f4 = Symbol('f4')
f5 = Symbol('f5')
f6 = Symbol('f6')
f7 = Symbol('f7')
f8 = Symbol('f8')
rho = Symbol('rho')
u = Symbol('u')
v = Symbol('v')

eq1 = f1+f5+f8-f3-f6-f7-rho*u
eq2 = f2+f5+f6-f4-f7-f8-rho*v
eq3 = f4-f2+2/3*rho*v

ans = solve([eq1,eq2,eq3],[f4,f7,f8])

pprint(ans)

# import numpy as np
# import matplotlib.pyplot as plt


# a = np.array([[3,3,3],[2,2,2],[1,1,1],[0,0,0]])
# # plt.imshow(a)
# # plt.colorbar()
# # plt.show()
# print(a[1:,1:])