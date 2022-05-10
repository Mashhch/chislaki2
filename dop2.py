import math
import matplotlib.pyplot as plt
import numpy as np

def next_2_ord(prev, tau,D, a, h):
    next = np.zeros(len(prev))
    for i in range(1, len(prev)-1):
        next[i] = - a*tau/(2*h)*(prev[i+1]-prev[i-1]) +1/2*(prev[i+1]-prev[i-1]) + D*tau/h**2*(prev[i+1]+ 2*prev[i] -prev[i-1])
    next[-1] = 0
    next[0] = 0
    return next

h = 0.25
a = 10
D = 0.5
N = 30
tau = min(h/a, h**2/(2*D))/2
x = np.arange(0, 10+h, h)
t = np.arange(0, N*tau,tau)
prev = [0]*len(x)
next = [0]*len(x)
for k in range(len(prev)):
    prev[k]= 0 if x[k]<1 or x[k]>2 else 1

for i in range(N):
    next = next_2_ord(prev, tau,D, a, h)
    prev = np.copy(next)

print(h/a, " ", h**2/(2*D), " ", tau)

plt.figure(3)
plt.plot(x, next, color='brown', label='u(x)')
plt.legend()
plt.grid(True)
plt.show()


