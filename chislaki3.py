import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


n = 1
k = -1/3
m = -1/3
alpha = 6**0.5
u0 = 1.25
hi = 1


def f(x):
    n = 1
    alpha = 6 ** 0.5
    if 0<= x and x<= alpha:
        return 1 - n/(2*n+4)*x**2
    else:
        return 0

def u_tochn(x, t):
    n = 1
    k = -1 / 3
    alpha = 6 ** 0.5
    u0 = 1.25
    hi = 1
    if t == 0:
        return 0
    elif 0<=x and x < alpha*(hi * u0**n * t**(k*n+1))**0.5:
        return u0 * t**k * (f(x/((hi * u0**n * t**(k*n+1))**0.5)))**(1/n)
    else:
        return 0

def lymbda(u):
    n = 1
    return u**n


def gamma1(t):
    k = -1 / 3
    u0 = 1.25
    return u0 * t**k

# n = 1
# k = 1
# m = 1
# alpha = 1
# u0 = 1.25
# hi = 1

# def f(x):
#     n = 1
#     alpha = 6 ** 0.5
#     if 0<= x and x<= alpha:
#         return 1 - x
#     else:
#         return 0
#
# def u_tochn(x, t):
#     n = 1
#     k = 1
#     alpha = 1
#     u0 = 1.25
#     hi = 1
#     if t == 0:
#         return 0
#     elif 0<=x and x <= alpha*(hi * u0**n * t**(k*n+1))**0.5:
#         return u0 * t**k * (f(x/((hi * u0**n * t**(k*n+1))**0.5)))**(1/n)
#     else:
#         return 0
#
# def lymbda(u):
#     n = 1
#     return u**n
#
#
# def gamma1(t):
#     k = 1
#     u0 = 1.25
#     return u0 * t**k

def gamma2(t):
    return 0


def phi(x):
    return 0

def coeffs():
    alpha = [0, 0]
    beta = [1, 1]
    return alpha, beta

def method_progonki(a,b,c,d, n):
    A = np.zeros(n)
    B = np.zeros(n)
    y = np.zeros(n)

    A[0] = -c[0] / b[0]
    B[0] = d[0] / b[0]

    for i in range(1, n - 1):
        A[i] = -c[i] / (b[i] + a[i] * A[i - 1])
    A[-1] = 0
    for i in range(1, n):
        B[i] = (d[i] - a[i] * B[i - 1]) / (b[i] + a[i] * A[i - 1])

    y[-1] = B[-1]
    for i in reversed(range(n-1)):
        y[i] = B[i] + A[i] * y[i + 1]
    return y



def next(prev, prev_iter, tau, t, x, h):
    alpha, beta = coeffs()
    a_const = 1

    a = np.zeros(len(x))
    b = np.zeros(len(x))
    c = np.zeros(len(x))
    d = np.zeros(len(x))

    d[0] = gamma1(t*tau)
    d[-1] = u_tochn(x[-1],tau*t)

    b[0] = -(alpha[0]/h) + beta[0]
    c[0] = alpha[0]/h

    b[-1] = (alpha[1] / h) + beta[1]
    a[-1] = -alpha[1]/h

    for i in range(1, len(x)-1):
        a[i] = tau/(2*h**2) * (lymbda(prev_iter[i]) + lymbda(prev_iter[i - 1]))
    for i in range(1, len(x)-1):
        b[i] = -tau/(2*h**2) * (2*lymbda(prev_iter[i]) + lymbda(prev_iter[i - 1]) + lymbda(prev_iter[i + 1])) - 1
    for i in range(1, len(x)-1):
        c[i] = tau/(2*h**2) * (lymbda(prev_iter[i]) + lymbda(prev_iter[i + 1]))
    for i in range(1, len(x)-1):
        d[i] = -prev[i]

    return method_progonki(a,b,c,d,len(x))





h = 0.01
x_left = 0
x_right = 5
x = np.arange(x_left, x_right+h, h)
tau = 0.01
u = np.zeros((2, len(x)))
u[0] = np.zeros(len(x))
u[1] = np.zeros(len(x))
err = np.zeros(len(x))

Number_frames = 200


for j in range(60):
    u[1] = next(u[0], u[1], tau, 1, x, h)

fig = plt.figure(1)
ax1 = plt.axes(xlim=(x_left-1, x_right), ylim=(-1, 3))
line1, = ax1.plot([], [], lw=3)
line2, = ax1.plot([], [], lw=3)
line1.set_data([], [])
line2.set_data([], [])


def animate(i):
    y2 = np.zeros(len(x))
    if i in [0, 1]:
        y1 = u[i]
        for j in range(len(x)):
            y2[j] = u_tochn(x[j], i * tau)
    else:
        u[0] = np.copy(u[1])
        for j in range(100):
            u[1] = next(u[0], u[1], tau, i, x, h)
        y1 = u[1]
        for j in range(len(x)):
            y2[j] = u_tochn(x[j], i * tau)

    line1.set_data(x, y1)
    line2.set_data(x, y2)
    line1.set_color("pink")
    line2.set_color("blue")
    return line2, line1



anim1 = FuncAnimation(fig, animate,
                     frames= 200, interval=100, blit=True)

anim1.save('iter.gif',  writer='pillow')

