import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import pandas as pd

df=pd.read_table('evolution.txt', sep='\s+',header=None)
fp=pd.read_table('double_well.txt',sep='\s+',header=None)

nx, ny=df.shape
mx,my=fp.shape
x=np.linspace(0.7213,1.823,nx)
y=df.iloc[0:nx,0]
px=fp.iloc[0:mx,0]
py=fp.iloc[0:mx,1]
#print(y)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
line,  = ax.plot(x, y)
ax1=ax.twinx()
power  = ax1.plot(px,py,'r--')
# ax1.set_yticks([500000,1000000])
# ax1.set_yticklabels(['5e5','1e6'])
ax.set_ylabel(r'$\Phi(x)$')
ax1.set_ylabel(r'$V(x)$')
ax.set_xlabel('x')
fig.legend(labels=['wave function','power function'],loc='upper right')

def animate(i):
	line.set_ydata(df.iloc[0:nx,i+1])
	return line,

def init():
	line.set_ydata(df.iloc[0:nx,1])
	power = ax1.plot(px,py)
	return power,line,

ani = animation.FuncAnimation(fig=fig, func=animate, frames=ny-1, interval=5, blit=False, repeat=True, repeat_delay= 1000)
ani.save('wave.gif', writer='imagemagick', fps=30)

#plt.show()

