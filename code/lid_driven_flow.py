import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

Lx = 100
Ly = 100

dL = 1
dT = 1

NX = Lx // dL
NY = Ly // dL

w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
e = [[0,0], [1,0], [0,1], [-1,0], [0,-1], [1,1], [-1,1], [-1,-1], [1,-1]]

X = np.arange(NX+1)
Y = np.arange(NY+1)
x,y = np.meshgrid(X,Y)

# plt.plot(x,y)
# plt.show()

rho = np.zeros((NX+1,NY+1))
u = np.zeros((NX+1,NY+1,2))
f = np.zeros((NX+1,NY+1,9))
f_eq = np.zeros((NX+1,NY+1,9))

c = dL / dT 
c_s = c / 3**0.5
rho_0 = 1 # density_0
U = 0.1 # blank speed
Re = 1000
nu = U * Lx / Re
tau = nu / c_s**2 + 0.5 * dT # relaxation time

rho[:,:] = rho_0
u[0,:,0] = U

for i in range(9):
    f[:,:,i] = w[i] * rho[:,:] * (1 + (u[:,:,0]*e[i][0] + u[:,:,1]*e[i][1])/c_s**2 +(u[:,:,0]*e[i][0] + u[:,:,1]*e[i][1])**2/(2*c_s**4) - (u[:,:,0]*u[:,:,0] + u[:,:,1]*u[:,:,1])/(2*c_s**2))

# lbm
MaxIteration = 100
for ite in range(MaxIteration):
    # collision
    for k in range(9):
        f_eq = w[k] * rho[:,:] * (1 + (u[:,:,0]*e[k][0] + u[:,:,1]*e[k][1])/c_s**2 +(u[:,:,0]*e[k][0] + u[:,:,1]*e[k][1])**2/(2*c_s**4) - (u[:,:,0]*u[:,:,0] + u[:,:,1]*u[:,:,1])/(2*c_s**2))
        f[:,:,k] += (f_eq - f[:,:,k]) * dT/tau
    # print('iteration:%d\n'%(ite+1),f,end='\n------------\n')
    # streaming
    f[:,1:,1] = f[:,:NY,1] # left -> right
    f[:,:NY,3] = f[:,1:,3] # right -> left
    f[:NX,:,2] = f[1:,:,2] # bottom -> top
    f[1:,:,4] = f[:NX,:,4] # top -> bottom
    f[:NX,1:,5] = f[1:,:NY,5] # leftbottom -> righttop
    f[1:,:NY,7] = f[:NX,1:,7] # righttop -> leftbottom
    f[1:,1:,8] = f[:NX,:NY,8] # lefttop -> rightbottom
    f[:NX,:NY,6] = f[1:,1:,6] # rightbottom -> lefttop

    # boundary condition
    # left bounce back
    f[:,0,1] = f[:,0,3]
    f[:,0,5] = f[:,0,7]
    f[:,0,8] = f[:,0,6]

    # right bounce back
    f[:,NY,3] = f[:,NY,1]
    f[:,NY,6] = f[:,NY,8]
    f[:,NY,7] = f[:,NY,5]

    # bottom bounce back
    f[NX,:,2] = f[NX,:,4]
    f[NX,:,5] = f[NX,:,7]
    f[NX,:,6] = f[NX,:,8]

    # moving lid
    rho_temp = f[0,1:NY,0]+f[0,1:NY,1]+f[0,1:NY,3]+2*(f[0,1:NY,2]+f[0,1:NY,5]+f[0,1:NY,6])
    f[0,1:NY,4] = f[0,1:NY,2]
    f[0,1:NY,7] = 0.5*f[0,1:NY,1]-0.5*f[0,1:NY,3]+f[0,1:NY,5]-0.5*rho_temp*U
    f[0,1:NY,8] = 0.5*f[0,1:NY,3]-0.5*f[0,1:NY,1]+f[0,1:NY,6]+0.5*rho_temp*U

    rho = np.sum(f,axis=2)

    u = np.zeros((NX+1,NY+1,2))
    for k in range(9):
        u[:,:,0] += f[:,:,k]*e[k][0]
        u[:,:,1] += f[:,:,k]*e[k][1]

    u[:,:,0] /= rho
    u[:,:,1] /= rho

    u[0,1:NY,0] = U
    u[0,1:NY,1] = 0

u_norm = (u[:,:,0]**2 + u[:,:,1]**2)**0.5
plt.imshow(u_norm)
plt.colorbar()
plt.show()