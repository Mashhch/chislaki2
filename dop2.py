import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter



def next_2_ord(prev, tau,D, a, h):
    next = np.zeros(len(prev))
    for i in range(1, len(prev)-1):
        next[i] = - a*tau/(2*h)*(prev[i+1]-prev[i-1]) +1/2*(prev[i+1]-prev[i-1]) + D*tau/h**2*(prev[i+1]- 2*prev[i] +prev[i-1])
    next[-1] = 0
    next[0] = 0
    return next


h = 0.025
a = 10
D = 0.25
x = np.arange(0, 10+h, h)
Number_frames = 200
tau = min(h/a, h**2/(2*D))/2
prev = [0]*len(x)
next = [0]*len(x)

for k in range(len(prev)):
    prev[k]= 0 if x[k]<1 or x[k]>2 else 1

next = next_2_ord(prev, tau,D, a, h)

print(h/a, " ", h**2/(2*D), " ", tau)


fig = plt.figure()
ax = plt.axes(xlim=(0, 10), ylim=(-3, 3))
line1, = ax.plot([], [], lw=3)

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
        next = next_2_ord(prev, tau,D, a, h)
        y = next

    line1.set_data(x, y)
    line1.set_color("yellow")
    return line1,

anim1 = FuncAnimation(fig, animate, init_func=init,
                      frames=Number_frames, interval=80, blit=True)

anim1.save('laba2_dop.gif', writer='pillow')
