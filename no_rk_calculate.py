import math
import numpy as np

show = True

# Constants
g = 9.81
m = 1
R = 0.1
L = 2*R
Ig = 1/12 * m * L**2
I0 = Ig + m * R**2
ue = 0.2
uc = 0.15

# Simulation parameters
output_file = "temp.txt"
dt = 0.00001
theta0 = math.pi/2 - 0.01
dtheta0 = 0

# Start file
file = open(output_file, 'w')
file.write("# Constants\n")
file.write("g = %f\n" % g)
file.write("m = %f\n" % m)
file.write("R = %f\n" % R)
file.write("L = %f\n" % L)
file.write("Ig = %f\n" % Ig)
file.write("I0 = %f\n" % I0)
file.write("ue = %f\n" % ue)
file.write("uc = %f\n" % uc)

file.write("\n# Simulation parameters\n")
file.write("dt = %f\n" % dt)
file.write("dtheta0 = %f\n" % dtheta0)

# General functions
def integrate(dt, x, dx):
    return x+dx*dt

def rk4(cte, dt, X0, f):
    k1 = dt*f(*cte, X0)
    k2 = dt*f(*cte, X0+k1/2.)
    k3 = dt*f(*cte, X0+k2/2.)
    k4 = dt*f(X0+k3)
    X = X0 + (k1+2*k2+2*k3+k4)/6.
    return X

def y_calc(R, theta):
    return R*math.sin(theta)

def dy_calc(R, theta, dtheta):
    return R*dtheta*math.cos(theta)

def ddy_calc(R, theta, dtheta, ddtheta):
    return R*(ddtheta*math.cos(theta)-dtheta**2*math.sin(theta))

def pos_calc(R, theta, x, y):
    dx = R*math.cos(theta)
    dy = R*math.sin(theta)
    return x-dx, x+dx, y-dy, y+dy

# No slip functions
def ii_noslip(dtheta, ddtheta):
    if (dtheta<0 or (dtheta==0 and ddtheta<0)):
        return 1
    else:
        return -1

def ftheta_noslip(cte, theta, dtheta, ddtheta_noslip):
    return np.array[dtheta, ddtheta_noslip(*cte, theta)]

def ddtheta_noslip(g, m, R, I0, theta):
    return -m*g*R/I0 * math.cos(theta)

def check_noslip(g, R, ue, theta, dtheta, ddtheta, ii):
    return not( ii*ddtheta >= -ii*(dtheta**2*math.cos(theta)+ue*(g/R-dtheta**2*math.sin(theta))*ii)/
                            (math.sin(theta)+ue*math.cos(theta)*ii) )

def N_noslip(g, m, R, theta, dtheta, ddtheta):
    return m*(g+R*(ddtheta*math.cos(theta)-dtheta**2*math.sin(theta)))

def x_noslip(R, dtheta):
    return R*math.cos(theta)

def dx_noslip(R, theta, dtheta):
    return -R*dtheta*math.sin(theta)

def ddx_noslip(R, theta, dtheta, ddtheta):
    return -R*(ddtheta*math.sin(theta)+dtheta**2*math.cos(theta))

# Slip functions
def ftheta_slip(cte, theta, dtheta, ddtheta_slip):
    return np.array[dtheta, ddtheta_slip(*cte, theta)]

def ddtheta_slip(g, m, R, Ig, theta, dtheta, ddtheta, ii):
    return ( (dtheta**2*math.sin(theta)-g/R)*(math.cos(theta)-uc*math.sin(theta)*ii)/
                (math.cos(theta)*(math.cos(theta)-uc*math.sin(theta)*ii)+Ig/(m*R**2)) )

def N_slip(g, m, R, theta, dtheta, ddtheta):
    return m*(g+R*(ddtheta*math.cos(theta)-dtheta**2*math.sin(theta)))

def ddx_slip(uc, N, ii):
    return uc*N*ii

def check_slip(R, theta, dtheta, dx, ii):
    return ii*dx < -dtheta*R*math.sin(theta)

def ii_slip(R, theta, dtheta, dx):
    if dx <= -dtheta*R*math.sin(theta):
        return 1
    else:
        return -1

# Simulation
t = 0
theta = theta0
dtheta = dtheta0

# Calculate if slip (assuming no slip)
ddtheta = ddtheta_noslip(g, m, R, I0, theta)

N = N_noslip(g, m, R, theta, dtheta, ddtheta)

ddx = ddx_noslip(R, theta, dtheta, ddtheta)
dx = dx_noslip(R, theta, dtheta)
x = x_noslip(R, theta)

ddy = ddy_calc(R, theta, dtheta, ddtheta)
dy = dy_calc(R, theta, dtheta)
y = y_calc(R, theta)

ii = ii_noslip(dtheta, ddtheta)
slip = check_noslip(g, R, ue, theta, dtheta, ddtheta, ii)

x1, x2, y1, y2 = pos_calc(R, theta, x, y)

if show:
    print(t, slip, ii, N, theta, dtheta, ddtheta, x, dx, ddx, y, dy, ddy, x1, y1, x2, y2)

file.write("\n# Simulation data\n")
file.write("# t \t slip \t ii \t N \t theta \t dtheta \t ddtheta \t x \t dx \t ddx \t y \t dy \t ddy \t x1 \t y1 \t x2 \t y2\n")
file.write("%f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" %
            (t, slip, ii, N, theta, dtheta, ddtheta, x, dx, ddx, y, dy, ddy, x1, y1, x2, y2))

while y>0:

    if not slip:
        ii = ii_noslip(dtheta, ddtheta)

        ddtheta = ddtheta_noslip(g, m, R, I0, theta)
        dtheta = integrate(dt, dtheta, ddtheta)
        theta = integrate(dt, theta, dtheta)

        N = N_noslip(g, m, R, theta, dtheta, ddtheta)

        ddx = ddx_noslip(R, theta, dtheta, ddtheta)
        dx = dx_noslip(R, theta, dtheta)
        x = integrate(dt, x, dx)

        ddy = ddy_calc(R, theta, dtheta, ddtheta)
        dy = dy_calc(R, theta, dtheta)
        y = y_calc(R, theta)

        slip = check_noslip(g, R, ue, theta, dtheta, ddtheta, ii)

    else:
        ii = ii_slip(R, theta, dtheta, dx)

        ddtheta = ddtheta_slip(g, m, R, Ig, theta, dtheta, ddtheta, ii)
        dtheta = integrate(dt, dtheta, ddtheta)
        theta = integrate(dt, theta, dtheta)

        N = N_slip(g, m, R, theta, dtheta, ddtheta)

        ddx = ddx_slip(uc, N, ii)
        dx = integrate(dt, dx, ddx)
        x = integrate(dt, x, dx)

        ddy = ddy_calc(R, theta, dtheta, ddtheta)
        dy = dy_calc(R, theta, dtheta)
        y = y_calc(R, theta)

        slip = check_slip(R, theta, dtheta, dx, ii)

    x1, x2, y1, y2 = pos_calc(R, theta, x, y)
    t += dt
    if show:
        print(t, slip, ii, N, theta, dtheta, ddtheta, x, dx, ddx, y, dy, ddy, x1, y1, x2, y2)

    file.write("%f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" %
                (t, slip, ii, N, theta, dtheta, ddtheta, x, dx, ddx, y, dy, ddy, x1, y1, x2, y2))
    #input()

file.close()
print("end")
