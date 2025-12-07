import numpy as np 
import os
import matplotlib.pyplot as plt
from src.SNR_utils import *
from src.bubbles   import SNR, SB
from src.utils 	   import set_axes_equal, colori
from src.shapes   import drawEllipsoid, draw_mysphere
import src.density_distributions as densities

Barros = densities.Barros()

# output dir 
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# initial conditions 

n   = 0.7       # cm^-3
xi   = 1.0      # metallicity in solar units

    
objects = []            
# ANTILLA
#objects += [bubble((0.0, 0.0, 0.0  ), 360.0, name="LB")]

# PER-TAU
objects += [SB((161.  ,-22.7, 218.0), 39.5, DensityProfile=Barros, name  = "PER-TAU",deltat_SNe=1e6)]
## GUM NEBULA
objects += [SB((258.  ,-2., 400.0), 36.0, DensityProfile=Barros, name   = "Gum Nebula")]
# GSH238_00_09
#objects += [SB((238.  ,-0., 790.0), 32.0, DensityProfile=Barros, name   = "GSH238_00_09")]
# Orion-Eridanus
objects += [SB((205.  ,-20., 290.0), 40.0, DensityProfile=Barros, name   = "Orion-Eridanus")]
objects += [SNR((275.5,18.4, 250  ), 23.0, DensityProfile=Barros, name   = "Antilla")]
# C1
objects += [SNR((7.   ,-18., 190.0), 15.0, DensityProfile=Barros, name   = "C1")]
# C2
objects += [SNR((12.  ,-11., 160.0), 33.0, DensityProfile=Barros, name   = "C2")]
## CHEPHEUS FLARE
objects += [SNR((120.  ,17., 300.0), 19.0, DensityProfile=Barros, name   = "Chepheus Flare")]

# MONOGEM RING 
objects += [SNR((203.  ,12., 300.0), 25.0, DensityProfile=Barros, name   = "Monogem Ring")]

# Vela - SNR
objects += [SNR((264.  ,-3.4, 290.0), 8.0, DensityProfile=Barros, name   = "Vela-SNR")]
# G354-33
trand =np.random.uniform()
objects += [SNR((354.  ,-33.5, 350.0), 14.0, DensityProfile=Barros, name = "G354-33",known_age=1e5*trand) ]
# G249+24
trand =np.random.uniform()
objects += [SNR((249.7 ,24.7 , 300.0), 4.5, DensityProfile=Barros, name  = "G249+24",known_age=1.5e4*trand) ]
# High latitude cavity 
objects += [SNR((150.0 ,50. , 307.0), 66.0, DensityProfile=Barros, name  = "h-l cavity") ]
# 2SECOND
objects += [SNR((148.27 ,33.72 , 207.0), 11.0, DensityProfile=Barros, name  = "h-l cavity2") ]

# NPCL
objects += [SNR((132.31 ,35.31 , 244.0), 25.0, DensityProfile=Barros, name  = "NPCL - NS #11") ]
objects += [SNR((133.76 ,32.11 , 488.0), 20.0, DensityProfile=Barros, name  = "NPCL - NS #12") ]
objects += [SNR((141.79 ,24.79 , 909.0), 10.0, DensityProfile=Barros, name  = "NPCL - NS #13") ]
# NPCL ALTERNATIVE

#objects += [SB((133.  ,32.11, 488.0), 60.0, name   = "npcÇl")]




# plots 
fig=plt.figure(0)

ax=fig.add_subplot(111, projection="3d")
colors=colori(len(objects),100)
eptimes=[]
f=open(file=f"{OUTPUT_DIR}/new_ic",mode="a")
f.write("# x,y,z, vx, vy, vz, mass, t_birth, dummy, particle_type, particle_tag, explosion_time\n")
for idx,obj in enumerate(objects):
    r     = obj.data_radius()
    t_age = obj.age(tmin=1e2, tmax=5e7, npoints=10000)/1e6
    obj.local_density(Barros)
    obj.Nsn_analytic()
    if obj.is_SB(): obj.add_random_SNe(radius=10)
    x,y,z = obj.helio()
    print(obj)
    tsn = 10
    tbirth=20.0-(obj.age()/1e6+tsn)
    texpl = tbirth+tsn
    if not obj.is_SB():
       f.write("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %d %d %d %10.5f\n"%(x,y,z,0.,0.,0.,21.6886, tbirth, 0, 6, -69,texpl))
    #    ax.scatter(x,y,z, color=colors[idx], s=10)
       eptimes+=[texpl]
    else:
       for snr in obj.progenitors:
           xi,yi,zi = snr.coord()
           tbirth=20.0-(snr.age()/1e6+tsn)
           texpl = tbirth+tsn


           f.write("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %d %d %d %10.5f\n"%(xi,yi,zi,0.,0.,0.,21.6886, tbirth, 0, 6, -69,texpl))
           ax.scatter(xi,yi,zi, color=colors[idx], s=10)
           eptimes+=[texpl]
    ax.text(x,y+10,z-20, obj.name, fontdict={"fontsize":8, "color":"pink"}, zorder=1000)
    #(xs,ys,zs) = drawSphere(x,y,z,r)
    
    ax = draw_mysphere(ax, x,y,z,r, color=colors[idx], ngrid= 100)
    #ax.plot_wireframe(xs, ys, zs, label=obj.name, lw=0, color=colors[idx])
    set_axes_equal(ax)


(xs,ys,zs) = drawEllipsoid((0,0,0),(150,150,250))
ax.plot_wireframe(xs, ys, zs, label="Local Bubble", color="red", rstride=10, alpha = 0.5)

(xs,ys,zs) = drawEllipsoid((300,0,0),(150,150,250))

#ax.plot_wireframe(xs, ys, zs, label="Loop I", color="blue", rstride=10)

# extra_points = np.array([
#     [-171.02956,   45.35633,  -13.22887],
#     [-171.02237,   67.44772,  -52.64560],
#     [ -87.61493, -414.26899,   -1.81436],
#     [ -65.18620, -446.10626,  -47.45587],
#     [ -42.42802, -410.61710,  -52.87772],
#     [-102.94746, -412.25627,  -59.63228],
#     [ -53.30997, -397.03501,  -21.10119],
#     [-176.50160, -124.84512, -122.74318],
#     [-230.13652, -152.87329, -119.69233],
#     [-318.72451, -118.50647,  -95.55306],
#     [  22.73643, -236.12689,   78.91226],
#     [ 179.35382,   22.02188,  -58.71323],
#     [ 153.62820,   32.65468,  -30.52944],
#     [-143.44571,  248.45526,   87.71151],
#     [-270.11688, -114.65781,   62.37351],
#     [ -30.25990, -287.90370,  -17.19885],
#     [ 290.26120,  -30.50768, -193.17794],
#     [ -94.55816, -255.62393,  125.36012],
# ])
# ax.scatter(extra_points[:,0], extra_points[:,1], extra_points[:,2], color='k', s=20, marker='x', label='SNe')


ax.set_xlabel("x [pc]")
ax.set_ylabel("y [pc]")
ax.set_zlabel("z [pc]")

ax.set_xlim(-500,500)
ax.set_ylim(-500,500)
ax.set_zlim(-500,500)

ax.legend(fontsize=8)
fig.savefig(f"{OUTPUT_DIR}/LB_neigh", dpi=512)
plt.show()
# fig=plt.figure(1)
# ax1=plt.subplot()
# ax1.scatter(eptimes,[0]*len(eptimes))
# plt.show()
