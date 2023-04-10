import matplotlib.pyplot as plt
import numpy as np


#alpha-Zr
pi = np.pi
deg = pi / 180

c = 5.149
a = 3.231

# dfpln = np.array([[1, 0, 1],
#                  [1, 0, 1],
#                  [1, 1, 0],
#                  [1, 1, 0]])  # дополнить всеми системами скольления данного семейства

#система скольжения {10-11}<11-20>
dfpln = np.array([[1, 0, -1, 1],
                  [1, -1, 0, 1],
                  [0, 1, -1, 1],
                  [-1, 1, 0, 1],
                  [-1, 0, 1, 1],
                  [0, -1, 1, 1]])

# dfdirct = np.array([[-1, 1, 1],
#                  [-1, -1, 1],
#                  [1, -1, 1],
#                  [-1, 1, 1]])    # дополнить всеми направлениями скольжения данного семейства

dfdirct = np.array([[1, -2, 1, 0],
                  [-1, -1, 2, 0],
                  [-2, 1, 1, 0],
                  [-1, -1, 2, 0],
                  [1, -2, 1, 0],
                  [2, -1, -1, 0]])

np.sum(dfpln * dfdirct, axis=1) 
                                    # в сумла квадратов по строкам, для проверки ортогональностик
                                      #    плоскости/направления§ внешняя система координат

X = np.array([np.sqrt(3) / 2, -0.5, 0])
Y = np.array([0, 1, 0])
Z = np.array([0, 0, 1]) # Определение направляющих косинусов плоскостей и направлений скольжений
                                        # во внешней системе координат по скалярному произведению

ar = [X, Y, Z]

pln_change_dim = np.array([[2 / np.sqrt(3), 1 / np.sqrt(3), 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 0, a / c]])

dirct_change_dim = np.array([[np.sqrt(3), np.sqrt(3) / 2, 0, 0],
                         [0, 1.5, 0, 0],
                         [0, 0, 0, c / a]])

dfpln = dfpln.T
dfdirct = dfdirct.T
dfpln = pln_change_dim.dot(dfpln)
dfdirct = dirct_change_dim.dot(dfdirct)
dfpln = dfpln.T
dfdirct = dfdirct.T


pln = np.array([]).reshape(0, 6)
for i in range(3):
  x = np.sum(dfpln*ar[i], axis=1)/np.sqrt(np.sum(dfpln**2,axis=1))
  pln = np.vstack([pln, x])

pln = pln.T
# pln

dirct = np.array([]).reshape(0, 6)
for i in range(3):
  x = np.sum(dfdirct*ar[i], axis=1)/np.sqrt(np.sum(dfdirct**2,axis=1))
  dirct = np.vstack([dirct, x])
dirct = dirct.T
# dirct

step = 5
nn = int(90 / step) + 1
phi1 = np.linspace(0, 0.5*pi, int(nn))
psi1 = np.linspace(0, 0.5*pi, int(nn))
psi, phi = np.meshgrid(psi1, phi1)

xx = np.tan(0.5*psi) * np.cos(phi)
yy = np.tan(0.5*psi) * np.sin(phi)
tdx = np.sin(psi) * np.cos(phi)
tdy = np.sin(psi) * np.sin(phi)
tdz = np.cos(psi)

# max(dfpln.shape)

tau = np.zeros((nn, nn))
for i in range(max(dfpln.shape)):
  chi = tdx*pln[i, 0] + tdy*pln[i, 1] + tdz*pln[i, 2]
  chi = chi / (np.sqrt(tdx**2 + tdy**2 + tdz**2) * np.sqrt(np.sum(pln[i,:]**2)))
  lambd = tdx*dirct[i, 0] + tdy*dirct[i, 1] + tdz*dirct[i, 2]
  lambd = lambd / (np.sqrt(tdx**2 + tdy**2 + tdz**2) * np.sqrt(np.sum(dirct[i,:]**2)))
  taul = abs(chi*lambd)
  indices = taul > tau
  tau[indices] = taul[indices]

# np.amax(tau)

# np.savetxt(sys.stdout, tau, delimiter='\t')

plt.contourf(xx, yy, tau, 11, cmap='coolwarm')
plt.colorbar()
plt.title("Фактор Шмида")
# plt.savefig('Факор Шмида 5 вариант.jpeg', dpi=300, format='jpeg')
plt.show()