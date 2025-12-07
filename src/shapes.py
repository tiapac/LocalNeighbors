import numpy as np 


def drawEllipsoid(center, radii, ngrid=20):
    u = np.linspace(0, 2 * np.pi, ngrid)
    v = np.linspace(0, np.pi, ngrid)
    x = center[0] + radii[0] * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radii[1] * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radii[2] * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z
    
def draw_mysphere(ax,xc,yc,zc,r,color, ngrid = 200):
    u = np.linspace(0, 2 * np.pi, ngrid)
    v = np.linspace(0, np.pi, ngrid)
    x = xc+r * np.outer(np.cos(u), np.sin(v))
    y = yc+r * np.outer(np.sin(u), np.sin(v))
    z = zc+r * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=0.5)#'linen'
    # plot circular curves over the surface
    theta = np.linspace(0, 2 * np.pi, ngrid)
    z = np.zeros(ngrid)+zc
    x = xc+r * np.sin(theta)
    y = yc+r * np.cos(theta)
    #z = np.zeros(100)+xc
    #y = zc+r * np.cos(theta)
    ax.plot(x, y, z, color='black', alpha=0.75)
    #ax.plot(z, x, y, color='black', alpha=0.75)

    ### add axis lines
    #zeros = np.zeros(1000)
    #line = np.linspace(-r,r,1000)

    #ax.plot(line, zeros, zeros, color='black', alpha=0.75)
    #ax.plot(zeros, line, zeros, color='black', alpha=0.75)
    #ax.plot(zeros, zeros, line, color='black', alpha=0.75)
    return ax
# Plot the surface