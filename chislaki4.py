import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def a(u, a_const, flag):
    if flag ==1:
        return a_const
    if flag ==2:
        return u


def f(u, a_const, flag):
    if flag ==1: return a_const*u
    if flag ==2: return (u**2)/2


def fi(x, flag):
    if flag ==1:
        return np.heaviside(x - 1, 1) * np.heaviside(2 - x, 1)
    if flag == 2:
        return (np.sin(np.pi * x / 2)) ** 2 * np.heaviside(x, 1) * np.heaviside(2 - x, 1)

def next_frame(prev, x,  h, tau, a_const, flag):

    next = np.zeros(len(x))
    for i in range(1, len(x) - 1):
        next[i] = 0.5*(prev[i+1] + prev[i-1]) - tau/(2*h)*(f(prev[i+1],a_const,flag) - f(prev[i-1], a_const,flag))
    next[0] = 0
    next[-1] = 0
    return next


a_const = 5
h = 0.01
C = 0.95
# flag = 1
flag =2
x = np.arange(0, 10+h, h)
Number_frames = 200

fig = plt.figure()
ax = plt.axes(xlim=(0, 10), ylim=(-3, 3))
line1, = ax.plot([], [], lw=3)



if flag == 1:
    prev = fi(x, flag)
    tau = (C * h) / a_const
    next = next_frame(prev, x, h, tau, a_const, 1)


else:
    prev = fi(x, 2)
    tau = (C * h) /max(a(prev, a_const, 2))
    next = next_frame(prev, x, h, tau,  a_const,2)


def init():
    line1.set_data([], [])
    return line1,


def animate(i):
    global prev,next
    if i == 0:
        y = prev
    elif i == 1:
        y = next
    else:
        prev = np.copy(next)
        # tau = (C * h) / a_const
        tau = (C * h) / max(a(prev, a_const, 2))
        # next = next_frame(prev, x, h, tau,  a_const, 1)
        next = next_frame(prev, x, h, tau, a_const, 2)
        y = next

    line1.set_data(x, y)
    line1.set_color("yellow")
    return line1,

anim1 = FuncAnimation(fig, animate, init_func=init,
                      frames=Number_frames, interval=50, blit=True)

anim1.save('number2.gif', writer='pillow')