import matplotlib.pyplot as plt
import numpy as np
# from pprint import pprint


def create_E_matrix(s, y1, y2, y3):
    E = (y1**4)*s[0,0] + (y2**4)*s[1,1] + (y3**4)*s[2,2]+ \
        + (y2**2)*(y3**2)*s[3,3] + y1**2*y3**2*s[4,4] + y1**2*y2**2*s[5,5]+ \
        + 2*s[0,1]*y1**2*y2**2 + 2*s[0,2]*y1**2*y3**2 + 2*s[1,2]*y2**2*y3**2
    E=1/E
    return E


def into_dec_cord(PSI1, PHI1, E_mat):
    X = E_mat * np.cos(PSI1) * np.cos(PHI1)
    Y = E_mat * np.cos(PSI1) * np.sin(PHI1)
    Z = E_mat * np.sin(PSI1)
    return X,Y,Z


def create_rot_matrix(ang1, ang2, ang3):
    rot_matrix = np.array([
                [np.cos(ang1)*np.cos(ang3)-np.sin(ang1)*np.sin(ang3)*np.cos(ang2),
                 np.cos(ang3)*np.sin(ang1)+np.sin(ang3)*np.cos(ang1)*np.cos(ang2),
                 np.sin(ang3)*np.sin(ang2)],
                [-np.cos(ang1)*np.sin(ang3)-np.sin(ang1)*np.cos(ang3)*np.cos(ang2),
                 -np.sin(ang1)*np.sin(ang3)+np.cos(ang1)*np.cos(ang3)*np.cos(ang2),
                 np.cos(ang3)*np.sin(ang2)],
                [np.sin(ang1)*np.sin(ang2), -np.cos(ang1)*np.sin(ang2), np.cos(ang2)]])
    return rot_matrix


def create_g4_matrix(cos_matrix):
    cos_matrix = np.array(cos_matrix)
    g4_matrix = np.array([
                [cos_matrix[0,0]**2, cos_matrix[1, 0]**2, cos_matrix[2,0]**2, 2*cos_matrix[1,0]*cos_matrix[2,0], 2*cos_matrix[0,0]*cos_matrix[2,0], 2*cos_matrix[1,0]*cos_matrix[0,0]],
                [cos_matrix[0, 1]**2, cos_matrix[1,1]**2, cos_matrix[2,1]**2, 2*cos_matrix[1,1]*cos_matrix[2,1], 2*cos_matrix[0,1]*cos_matrix[2,1], 2*cos_matrix[1,1]*cos_matrix[0,1]],
                [cos_matrix[0, 2]**2, cos_matrix[1,2]**2, cos_matrix[2,2]**2, 2*cos_matrix[1,2]*cos_matrix[2,2], 2*cos_matrix[0,2]*cos_matrix[2,2], 2*cos_matrix[1,2]*cos_matrix[0,2]],
                [cos_matrix[0, 1]*cos_matrix[0,2], cos_matrix[1,1]*cos_matrix[1,2], cos_matrix[2,1]*cos_matrix[2,2], cos_matrix[1,1]*cos_matrix[2,2] + cos_matrix[1,2]*cos_matrix[2,1], cos_matrix[2,1]*cos_matrix[0,2] + cos_matrix[2,2]*cos_matrix[0,1], cos_matrix[1,1]*cos_matrix[0,2] + cos_matrix[1,2]*cos_matrix[0,1]],
                [cos_matrix[0, 2]*cos_matrix[0,0], cos_matrix[1,2]*cos_matrix[1,0], cos_matrix[2,2]*cos_matrix[2,0], cos_matrix[1,2]*cos_matrix[2,0] + cos_matrix[1,0]*cos_matrix[2,2], cos_matrix[2,2]*cos_matrix[0,0] + cos_matrix[2,0]*cos_matrix[0,2], cos_matrix[1,2]*cos_matrix[0,0] + cos_matrix[1,0]*cos_matrix[0,2]],
                [cos_matrix[0, 0]*cos_matrix[0,1], cos_matrix[1,0]*cos_matrix[1,1], cos_matrix[2,0]*cos_matrix[2,1], cos_matrix[1,0]*cos_matrix[2,1] + cos_matrix[1,1]*cos_matrix[2,0], cos_matrix[2,0]*cos_matrix[0,1] + cos_matrix[2,1]*cos_matrix[0,0], cos_matrix[1,0]*cos_matrix[0,1] + cos_matrix[1,1]*cos_matrix[0,0]]])
    return g4_matrix


def rotate_C(c, ang1, ang2, ang3):
    rot_mat = create_rot_matrix(ang1, ang2, ang3)
    guid_cos_matrix = [[rot_mat[i,j] / np.sqrt(rot_mat[0,j]**2 + rot_mat[1,j]**2 + rot_mat[2,j]**2) for j in range(3)]
                       for i in range(3)]
    g4 = create_g4_matrix(guid_cos_matrix)
    C1 = np.zeros((6, 6))
    for k in range(6):
        for l in range(6):
            for m in range(6):
                for n in range(6):
                    C1[k,l] += c[n,m]*g4[k,n]*g4[l,m]
    return C1

pi = np.pi
deg = np.pi/180
c11 = 350
c12 = 67
c44 = 101
gran1 = [0, 0, 0]
gran2 = [0, 0, 0]
gran3 = [45, 45, 0]
gran4 = [45, 55, 0]
gran5 = [60, 0, 0]

C = np.array([[c11, c12, c12, 0, 0, 0],
              [c12, c11, c12, 0, 0, 0],
              [c12, c12, c11, 0, 0 ,0],
              [0, 0, 0, c44, 0, 0],
              [0, 0, 0, 0, c44, 0],
              [0, 0, 0, 0, 0, c44]])

S = np.linalg.inv(C)
Grains = np.array([gran1,
                  gran2,
                  gran3,
                  gran4,
                  gran5]) * deg

phi = np.linspace(0,2*pi, 100)
psi = np.linspace(-pi,pi, 100)
PHI, PSI = np.meshgrid(phi, psi)

#Монокристалл

x = np.cos(PHI)*np.cos(PSI)
y = np.sin(PHI)*np.cos(PSI)
z = np.sin(PSI)
E1 = create_E_matrix(S, x, y, z)
dec_cord1 = into_dec_cord(PSI, PHI, E1)

#Поликристалл

X1 = 0
Y1 = 0
Z1 = 0
for i in range(0,5):
    C_rotated = rotate_C(C, Grains[i,0], Grains[i,1], Grains[i,2])
    S1 = np.linalg.inv(C_rotated)
    E2 = create_E_matrix(S1, x, y, z)
    dec_cord2 = into_dec_cord(PSI, PHI, E2)
    X1 = (X1*(i) + dec_cord2[0]) / (i + 1)
    Y1 = (Y1*(i) + dec_cord2[1]) / (i + 1)
    Z1 = (Z1*(i) + dec_cord2[2]) / (i + 1)

#Создаю картинку

fig = plt.figure()
fig.suptitle('Распределение модуля упругости Е для W')
ax = fig.add_subplot(121, projection='3d')
cmap = plt.get_cmap('inferno')
norm = plt.Normalize(dec_cord1[2].min(), dec_cord1[2].max())
facecolors = cmap(norm(dec_cord1[2]))
surf1 = ax.plot_surface(dec_cord1[0], dec_cord1[1], dec_cord1[2], facecolors=facecolors,
                        shade=False, rstride=1, cstride=1)
surf1.set_edgecolors('k')
surf1.set_linewidth(0.2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Монокристалл')
ax1 = fig.add_subplot(122, projection='3d')
norm1 = plt.Normalize(Z1.min(), Z1.max())
surf1.set_edgecolors('k')
surf1.set_linewidth(0.2)
facecolors = cmap(norm(Z1))
surf2 = ax1.plot_surface(X1, Y1, Z1, facecolors=facecolors,
                         shade=False, rstride=1, cstride=1)
surf2.set_edgecolors('k')
surf2.set_linewidth(0.2)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.set_title('Поликристалл 5 зерен')
plt.show()

# v = np.amax(E1, axis=0)
# v = v.reshape(-1,1)
# pprint(v)
