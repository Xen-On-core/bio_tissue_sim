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


def g(c):
    if c >= 1:
        return 1
    elif c >= 0:
        return c
    else:
        return 0


def div_a(x):
    if x <= 0:
        return 0
    elif x >= 1:
        return 1
    else:
        return x


time_scale = 10
T = 100 * time_scale
M = 100
h = 1 / M
tau = 1 / T * time_scale
A = 5
B = 5
D = 5

point = np.arange(0, 1 + h, h)

s = [[0 for _ in range(M)] for _ in range(T)]
phi = [[0 for _ in range(M)] for _ in range(T)]
phi0 = [0.5 for _ in range(M)]
C = [[0 for _ in range(M)] for _ in range(T)]
C0 = [0.5 for _ in range(M)]
for i in range(M):
    # s[0][i] = np.exp(-(point[i]+0.5)**10)/2
    # phi[0][i] = np.exp(-(point[i]-1.5)**4)/2
    s[0][i] = point[i] * (1 - point[i])
    phi[0][i] = 0.5
    C[0][i] = 0.5

alpha = [0 for _ in range(M + 1)]
beta = [0 for _ in range(M + 1)]

# for t in range(T - 1):
#     alpha[1] = 1
#     beta[1] = 0
#
#     for i in range(1, M - 1):
#         f1 = B * (phi[t][i + 1] * s[t][i + 1] + phi[t][i] * s[t][i]) / 2
#         f2 = B * (phi[t][i - 1] * s[t][i - 1] + phi[t][i] * s[t][i]) / 2
#         a = tau * f1 / (h ** 2 * phi[t][i])
#         b = tau * f2 / (h ** 2 * phi[t][i])
#         c = 1 + a + b
#         f = tau * F_s(phi[t][i], s[t][i]) + s[t][i]
#         alpha[i+1] = b / (c - alpha[i] * a)
#         beta[i+1] = (beta[i] * a + f) / (c - alpha[i] * a)
#
#     s[t + 1][M - 1] = beta[M - 1] / (1 - alpha[M - 1])
#     # s[t + 1][M - 1] = 0.5
#     for i in range(M - 2, -1, -1):
#         s[t + 1][i] = alpha[i + 1] * s[t + 1][i + 1] + beta[i + 1]
#
#     for i in range(M):
#         k1 = F(phi[t][i], s[t + 1][i])
#         k2 = F(phi[t][i] + tau / 2 * k1, s[t + 1][i])
#         k3 = F(phi[t][i] + tau / 2 * k2, s[t + 1][i])
#         k4 = F(phi[t][i] + tau * k3, s[t + 1][i])
#         phi[t + 1][i] = phi[t][i] + tau / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


for t in range(T - 1):
    alpha[1] = 0
    beta[1] = 0

    for i in range(1, M - 1):
        # f1 = B * (phi0[i + 1] * s[t][i + 1] + phi0[i] * s[t][i]) / 2
        f1 = B * (phi0[i + 1] * div_a(s[t][i + 1]) + phi0[i] * div_a(s[t][i])) / 2
        # f2 = B * (phi0[i - 1] * s[t][i - 1] + phi0[i] * s[t][i]) / 2
        f2 = B * (phi0[i - 1] * div_a(s[t][i - 1]) + phi0[i] * div_a(s[t][i])) / 2
        a = tau * f1 / (h ** 2)
        b = tau * f2 / (h ** 2)
        c = phi0[i] + a + b
        f = tau * Ft(phi0[i], s[t][i]) + phi0[i] * s[t][i]
        alpha[i + 1] = b / (c - alpha[i] * a)
        beta[i + 1] = (beta[i] * a + f) / (c - alpha[i] * a)

    s[t + 1][M - 1] = 0  # beta[M - 1] / (1 - alpha[M - 1])
    # s[t + 1][M - 1] = 0.5
    for i in range(M - 2, -1, -1):
        s[t + 1][i] = alpha[i + 1] * s[t + 1][i + 1] + beta[i + 1]

    alpha[1] = 0
    beta[1] = 0.5

    for i in range(1, M - 1):
        a = tau * D / (h ** 2)
        b = tau * D / (h ** 2)
        c = 1 + a + b
        f = C[t][i] - tau * Ft(phi0[i], s[t][i])*g(C[t][i])
        alpha[i + 1] = b / (c - alpha[i] * a)
        beta[i + 1] = (beta[i] * a + f) / (c - alpha[i] * a)

    C[t + 1][M - 1] = 0
    for i in range(M - 2, -1, -1):
        C[t + 1][i] = alpha[i + 1] * C[t + 1][i + 1] + beta[i + 1]


def view3d():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    X = np.linspace(0, 1, M)
    Y = np.linspace(0, 1 * time_scale, T)
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
    ax.plot(s[T - 1][:-1])
    ax.plot(phi[T - 1][:-1])
    print(np.max(np.array(s[T - 1][:-1])))
    print(np.max(np.array(phi[T - 1][:-1])))


view3d()

print('\n')

print(np.max(np.array(s)))
print(np.min(np.array(s)))
# print(np.max(np.array(phi)))
# print(np.min(np.array(phi)))

plt.show()
