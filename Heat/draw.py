import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import pandas as pd

df=pd.read_table('heat_1.txt', sep='\s+',header=None)

#print (df.iloc[0,0])
df2=pd.read_table('variance.txt', sep='\s+')
#exit()
mx, my=df2.shape
#boundary=df2.iloc[mx-1,1]
nx, ny=df.shape
x=np.arange(0,nx,1)
y=df.iloc[0:nx,0]
#print(y)
fig, ax = plt.subplots()
line,   = ax.plot(x, y)
ax.set_xlabel(r'x/mm')
ax.set_ylabel('temperature/'+r'$^\circ $C')
ax.set_xticks([0,20,40,60,80,100,120,140,160])
ax.set_xticklabels(['0','2','4','6','8','10','12','14','16'])

def animate(i):
    if i>mx:
        line.set_color('r')
    line.set_ydata(df.iloc[0:nx,i+1])
    return line,

def init():
    line.set_ydata(df.iloc[0:nx,1])
    return line,

ani = animation.FuncAnimation(fig=fig, func=animate, frames=ny-1, interval=2, blit=False, repeat=True, repeat_delay= 1000)
ani.save('heat.gif', writer='imagemagick', fps=30)

#plt.show()
