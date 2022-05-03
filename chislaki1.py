import math
import matplotlib.pyplot as plt
import numpy as np


def f(x, t):
    return 5 * (x ** 4 - 1 + x * t * (t + 2 * x) ** 2) / (3 * (1 + x ** 2 * (x + t) ** 2)**2)


def fi1(x):
    return 5 / 3 * np.arctan(x ** 2)

def fi1_xx(x):
    return (-10/3) *(4 * x**4/(x**4+1)- 1)/(x**4+1)

def fi2(x):
    return 5 * x / (3 * (1 + x ** 4))


def coefs():
    alpha = [-1, 1]
    beta = [2, 0]
    return alpha, beta


def gamma1(t):
    return -5 * t / 3


def gamma2(t):
    return 5 * (2 + t) / (3 * (1 + (1 + t) ** 2))


def u0(x, t):
    return 5 / 3 * np.arctan(x ** 2 + x * t)


def u_1_ord(x_left, x_right, h, t0, T, t):
    alpha, beta = coefs()
    x_knots = np.arange(x_left, x_right+h, h)
    t_knots = np.arange(t0, T, t)

    u = [[0] * len(x_knots) for i in range(3)]
    # первая координата время
    for j in range(len(x_knots)):
        u[0][j] = fi1(x_knots[j])
        u[1][j] = fi2(x_knots[j]) * t + fi1(x_knots[j])

    for i in range(2, len(t_knots)):
        for k in range(1, len(x_knots) - 1):
            u[2][k] = (a * t / h) ** 2 * (u[1][k + 1] - 2 * u[1][k] + u[1][k - 1]) + t ** 2 * f(x_knots[k], t_knots[i - 1]) + 2 * u[1][k] - u[0][k]

        u[2][0] = (gamma1(t_knots[i]) - (alpha[0] / h) * u[2][1]) * (h / (beta[0] * h - alpha[0]))
        u[2][-1] = (gamma2(t_knots[i]) + (alpha[1] / h) * u[2][-2]) * (h / (beta[1] * h + alpha[1]))

        u[0] = [u[1][k] for k in range(len(u[0]))]
        u[1] = [u[2][k] for k in range(len(u[0]))]

    z0 = [0] * len(x_knots)
    for k in range(len(x_knots)):
        z0[k] = u0(x_knots[k], t_knots[-1])

    return u[-1], z0

def u_2_ord(x_left, x_right, h, t0, T, t):
    alpha, beta = coefs()
    x_knots = np.arange(x_left, x_right+h, h)
    t_knots = np.arange(t0, T, t)
    z0 = [0] * len(x_knots)
    for k in range(len(x_knots)):
        z0[k] = u0(x_knots[k], t_knots[-1])
    u = [[0] * len(x_knots) for i in range(3)]
    # первая координата время
    u[0] = fi1(x_knots)
    u[1] = fi2(x_knots) * t + fi1(x_knots) + (t ** 2 / 2) * (0.5 * (-10/3) *(4 * x_knots**4/(x_knots**4+1)- 1)/(x_knots**4+1) + f(x_knots, 0))

    for i in range(2, len(t_knots)):
        for k in range(1, len(x_knots) - 1):
            u[2][k] = 0.5*( t / h) ** 2 * (u[1][k + 1] - 2 * u[1][k] + u[1][k - 1]) + t ** 2 * f(k*h, t*i-t) + 2 * u[1][k] - u[0][k]

        u[2][0] = 2*0.5*( t / h) ** 2 *(u[1][1] -(1 - h*beta[0]/alpha[0])*u[1][0] - h/alpha[0] * gamma1(i*t-t)) + 2*u[1][0] - u[0][0] + t**2*f(0, t*i-t)
        u[2][-1] = 2*0.5*( t / h) ** 2 *(u[1][-2] -(1 + h*beta[1]/alpha[1])*u[1][-1] + h/alpha[1] * gamma2(i*t-t)) + 2*u[1][-1] - u[0][-1] + t**2*f(x_knots[-1], t*i-t)

        u[0] = [u[1][k] for k in range(len(u[0]))]
        u[1] = [u[2][k] for k in range(len(u[0]))]

    return u[-1], z0

x_left = 0
x_right = 1
a = 1/np.sqrt(2)
t0 = 0
T = 1

#matrix
t = 0.005
h = 0.01
x_ = np.arange(x_left, x_right+h, h)
z_1, z0 = u_1_ord(x_left, x_right, h, t0, T, t)
y1 = [abs(z_1[k]-z0[k]) for k in range (len(z_1))]
er111 = max(y1)
print(er111)
z_2, z02 = u_2_ord(x_left, x_right, h, t0, T, t)
print("lol")
y2 = [abs(z_2[k]-z02[k]) for k in range (len(z_2))]
er222 = max(y2)
print(er222)

plt.figure(1)
plt.subplot(211)
plt.plot(x_, z0, color='yellow', label='Истинное значение')
plt.plot(x_, z_1, color='green', label='1 порядка')
plt.title(' график функции')
plt.legend()
plt.grid(True)

plt.subplot(212)
plt.plot(x_, y1, color='yellow', label='Ошибка 1 порядка')
plt.plot(x_, [h]*len(x_), color='green', label='h')
plt.legend()
plt.grid(True)

plt.figure(2)
plt.subplot(211)
plt.plot(x_, z0, color='yellow', label='Истинное значение')
plt.plot(x_, z_2, color='brown', label='2 порядка')
plt.title(' график функции')
plt.legend()
plt.grid(True)

plt.subplot(212)
plt.plot(x_, y2, color='yellow', label='Ошибка 2 порядка')
plt.plot(x_, [h**2]*len(x_), color='green', label='h^2')
plt.legend()
plt.grid(True)


# hrange = np.arange(hmin, hmax, hstep)
hrange = [0.05, 0.005, 0.025, 0.125, 0.1, 0.002]
h_log = np.zeros(len(hrange))
error_1 = np.zeros(len(hrange))
error_2 = np.zeros(len(hrange))
for j in range(len(hrange)):
    h = hrange[j]
    t = h/5
    t_ = np.arange(t0, T, t)
    h_log[j] = np.log(h)
    z_number_1, z_true = u_1_ord(x_left, x_right, h, t0, T, t)
    y1 = [abs(z_number_1[k] - z_true[k]) for k in range(len(z_true))]
    error_1[j] = np.log(max(y1))
    z_number_2, z_true_2 = u_2_ord(x_left, x_right, h, t0, T, t)
    y2 = [abs(z_number_2[k] - z_true_2[k]) for k in range(len(z_true_2))]
    error_2[j] = np.log(max(y2))

plt.figure(3)
plt.plot(h_log, error_1, color='yellow', label='1 порядка')
plt.plot(h_log, error_2, color='brown', label='2 порядка')
plt.title('Логарифмическая ошибка')
plt.legend()
plt.grid(True)
plt.show()