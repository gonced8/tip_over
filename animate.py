import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import patches

interactive=True
gif=False
video=False

speed=1./8
FPS=60

filename = "temp.txt"

# Reading constants
with open(filename, 'r') as file:
    line = file.readline()
    cte = []
    for i in range(8):
        line = file.readline()
        line = line.split()
        cte += [float(line[-1])]
    cte = tuple(cte)

# Reading simulation data
data = np.loadtxt(filename, skiprows=16, usecols=(0, 1, 4, 7, 10, 13, 14, 15, 16))

data[-1, :] = data[-2, :] #ignore last point because theta<0
t = data[:, 0]
slip = data[:, 1].astype(bool)
theta = data[:, 2]
x = data[:, 3]
y = data[:, 4]
p1 = data[:, 5:7]
p2 = data[:, 7:9]

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(9, 5)

L = 2*cte[2]
ax = plt.axes(xlim=(-L*1.1, L*1.1), ylim=(-L*0.1, L*1.1))

ground = patches.FancyArrowPatch(posA=(-L*1.1, 0), posB=(L*1.1, 0), ls='solid', lw=3., color='black')
body = patches.FancyArrowPatch(posA=p1[-1], posB=p2[-1], ls='solid', lw=5.)
trail, = ax.plot([], [], lw=1, c='red')

if interactive:

    def init():
        if slip[0]:
            body.set_color('g')
        else:
            body.set_color('b')
        body.set_positions(p1[0], p2[0])
        ax.add_patch(ground)
        ax.add_patch(body)
        trail.set_data(x[0], y[0])
        return ground, body, trail,

    def animate(i, step, frames):
        if (i+1)==frames:
            i=-2
            step=1

        if slip[i*step]:
            body.set_color('g')
        else:
            body.set_color('b')
        body.set_positions(p1[i*step], p2[i*step])
        trail.set_data(x[:i*step+1], y[:i*step+1])
        return ground, body, trail,

    interval=int(1000./FPS)
    if interval==0:
        print("Interval set to 1ms")
        interval=1
    step=int(speed/((t[1]-t[0])*FPS))
    frames=int((len(t))/step+0.5)

    anim = animation.FuncAnimation(fig, animate, fargs=(step, frames,),
                                   init_func=init,
                                   frames=frames,
                                   interval=interval,
                                   blit=True, repeat=True)

else:
    ax.add_patch(ground)
    ax.add_patch(body)
    trail.set_data(data[:,-2], data[:,-1])

if gif:
    anim.save('./gif.gif', writer='imagemagick', fps=FPS)

if video:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=FPS, metadata=dict(artist='Goncalo Raposo'), bitrate=5000)
    anim.save('./mp4.mp4', writer=writer)

ax.set_aspect('equal')
plt.tight_layout()
plt.show()

print('end')
