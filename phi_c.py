import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib


def Ft(x, y):
    return A * y * (1 - y) * x * (1 - x)


def F(x, y):
    return A * (1 - y) * x * (1 - x)


def F_s(x, y):
    return A * (1 - y) ** 2 * (1 - x)


def g(x):
    if x >= 1:
        return 1
    elif x >= 0:
        return x
    else:
        return 0


def div_a(x):
    if x <= 0:
        return 0
    elif x >= 1:
        return 1
    else:
        return x


time_scale = 2
space_scale = 5
T = 100 * time_scale
M = 100 * space_scale
h = 1 / M * space_scale
tau = 1 / T * time_scale
A = 1
B = 1
D = 1
hi = 1

point = np.arange(0, space_scale + h, h)

s = [[0 for _ in range(M)] for _ in range(T)]
phi = [[0 for _ in range(M)] for _ in range(T)]
phi0 = [0.5 for _ in range(M)]
C = [[0 for _ in range(M)] for _ in range(T)]
C0 = [0.5 for _ in range(M)]
for i in range(M):
    s[0][i] = 0.5
    phi[0][i] = 0.5
    C[0][i] = 0.9/(i+1)

alpha = [0 for _ in range(M + 1)]
beta = [0 for _ in range(M + 1)]

for t in range(T - 1):
    alpha[1] = 0
    beta[1] = 0.9

    for i in range(1, M):
        a = tau * D / (h ** 2)
        b = tau * D / (h ** 2)
        c = 1 + a + b
        f = C[t][i] - tau*A * Ft(phi0[i], s[t][i])*g(C[t][i])
        alpha[i + 1] = b / (c - alpha[i] * a)
        beta[i + 1] = (beta[i] * a + f) / (c - alpha[i] * a)

    C[t+1][M-1] = 0.9/(M)
    for i in range(M - 2, -1, -1):
        C[t + 1][i] = alpha[i + 1] * C[t + 1][i + 1] + beta[i + 1]

    alpha[1] = 1
    beta[1] = 0

    for i in range(1, M - 1):
        f1 = B * (phi0[i + 1] * div_a(s[t][i + 1]) + phi0[i] * div_a(s[t][i])) / 2
        f2 = B * (phi0[i - 1] * div_a(s[t][i - 1]) + phi0[i] * div_a(s[t][i])) / 2
        a = tau * f1 / (h ** 2)
        b = tau * f2 / (h ** 2)
        c = phi0[i] + a + b
        f = tau * (F_s(phi[t][i], s[t][i])*g(C[t][i]) + hi*f1*(C[t+1][i+1] - C[t+1][i]) - hi*f2*(C[t+1][i] - C[t+1][i-1])) + s[t][i]*phi0[i]
        alpha[i+1] = b / (c - alpha[i] * a)
        beta[i+1] = (beta[i] * a + f) / (c - alpha[i] * a)

    s[t + 1][M - 1] = beta[M - 1] / (1 - alpha[M - 1])
    # s[t + 1][M - 1] = 0.5
    for i in range(M - 2, -1, -1):
        s[t + 1][i] = alpha[i + 1] * s[t + 1][i + 1] + beta[i + 1]


def view3d():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X = np.linspace(0, space_scale, M)
    Y = np.linspace(0, time_scale, T)
    X, Y = np.meshgrid(X, Y)

    ax.set_xlabel('$x$')
    ax.set_ylabel('$t$')

    surf1 = ax.plot_surface(X, Y, np.array(s), label='$s$')
    surf1._facecolors2d = surf1._facecolor3d
    surf1._edgecolors2d = surf1._edgecolor3d

    # surf2 = ax.plot_surface(X, Y, np.array(phi), label='$\phi$')
    # surf2._facecolors2d = surf2._facecolor3d
    # surf2._edgecolors2d = surf2._edgecolor3d

    surf3 = ax.plot_surface(X, Y, np.array(C), label='$C$')
    surf3._facecolors2d = surf3._facecolor3d
    surf3._edgecolors2d = surf3._edgecolor3d

    ax.legend()


def view2d():
    fig, ax = plt.subplots()
    ax.plot(point[:-2], s[T - 1][:-1])
    ax.plot(point[:-2], C[T - 1][:-1])


view3d()

plt.show()
