import itasca as it, math as ma, numpy as np, bisect as bi, random as ra, time as ti
from vec import vec

it.command("python-reset-state false")

# define platter-mesh contact parameters
kn = 1.0e9
fric = 0.2
ks = kn

# pullout velocity
wallVel = 10

# mesh size
net_dim = 3.0

# apostrophe
ap = "'"

# Numerical parameters from Thoeni et al., 2013
net_ball_radius = 0.0108/2 # m
brad = net_ball_radius
net_ball_mass = 4.3*1e-3 # kg
net_ball_density = net_ball_mass/(4.0/3*ma.pi*net_ball_radius**3)

it.command("model new")    
it.command("program load contactmodelmechanical '../contactmodelmechanical3dmacroelement007.dll'")
it.command("model domain extent -10 10")

# parameters for the lengths of individual bits of mesh (dt = double twist, sw = single wire extent in the direction of the double-twist, mos = in the orthogonal direction)
mos = 0.08
dt = 0.04
sw =  0.03

sw_gap = (dt**2+sw**2)**0.5 - net_ball_radius*2
dt_gap = dt - net_ball_radius*2

# radius multiplier factors
rmuldt =  dt_radius/net_ball_radius
rmulsw =  single_wire_diameter/2/net_ball_radius


# Numerical parameters from Thoeni et al., 2013
net_ball_radius = 0.0108/2 # m
net_ball_mass = 4.3*1e-3 #[kg]
net_ball_density = net_ball_mass/(4.0/3*ma.pi*net_ball_radius**3)


b = sw
a = dt


# generate the mesh as a repeating pattern (script structure taken from the equivalent one in YADE)
net_starting = -net_dim/2
xstart = net_starting-mos
ystart = net_starting
z = 0
ab = a + b
ny = int( (net_dim-a)/ab ) + 1
ly= ny*a+(ny-1)*b
nx = int(net_dim/mos) + 1
lx = (nx-1)*mos
idx  = 0
for i in range(0,ny+1):
    y = ystart + i*ab
    for j in range(1,nx+1):
        x = xstart + j*mos
        ballpos = vec([x,y,z])
        strpos = ' '.join(str(s) for s in ballpos)
        it.command("ball create radius " + str(net_ball_radius) + " group "+ap+'dt'+str(idx)+ap+" position " + strpos)
        ballpos = vec([x,y+a,z])
        strpos = ' '.join(str(s) for s in ballpos)
        it.command("ball create radius " + str(net_ball_radius) + " group "+ap+'dt'+str(idx)+ap+" position " + strpos)
    xstart = xstart - 0.5 * mos*pow(-1,i+1)
    idx += 1
    idx = ma.ceil(idx%2)
    nx = int(nx+1*pow(-1,i+1))

# load the macroelement properties
it.command("program call 'load_properties'")

## Instantiate the contacts  
it.command("""
contact cmat proximity 0.05
contact cmat apply
model clean
contact property lin_mode 1
contact method bond gap 4
model energy mechanical on
model large-strain on
model orientation-tracking on
""")


# fix the mesh extremities (this is only done in one direction here, do the same in X if you want to
it.command("ball group 'lock=bottom' range position-y -99 "+str(xstart+dt*2+1e-3))
it.command("ball group 'lock=top' range position-y "+str(xstart+nx*mos-1e-3)+" 99")

it.command("ball fix velocity range group 'lock=bottom' union group 'lock=top'")
num_damping_net = 0.2
it.command("ball attribute density "+str(net_ball_density)+" damp "+str(num_damping_net))
it.command("model gravity 9.81")


# add ball-facet contacts for mesh-punch indenter contacts
it.command("contact cmat default type ball-facet model linear property kn "+str(kn)+" ks "+str(ks)+" fric "+str(fric)) #ks = 0.1 kn according to thoeni et al., 2013

# deactivate ball-ball contact detection (saves computational cost, we only care about mesh-platter contacts)
it.command("""
contact cmat remove all
contact detection off type ball-ball
""")

# import platter
it.command("wall import filename 'punching_test_geometry.stl'")

# place it just below the lowest particle
myWall = it.wall.find(1)
initialWallZ = myWall.pos()[2]
maxzw = -999
minzw = 999
for v in myWall.vertices():
    maxzw = max([maxzw,v.pos()[2]])
    minzw = min([minzw,v.pos()[2]])
minBz = 999
for b in it.ball.list():
    minBz = min([b.pos()[2],minBz])

verticalWallExtent = maxzw-minzw
verticalWallOffset = maxzw - initialWallZ

myWall.set_pos([0,0,-verticalWallOffset - net_ball_radius+ minBz])
myWall.set_vel((0,0,wallVel))

# calculate force-displacement variables
def forceDisp():
    if it.cycle()%output_frequency==0:
        Force = np.linalg.norm(myWall.force_contact())
        Disp = wallVel * (it.mech_age()-t0)
        it.fish.set('Force',Force)
        it.fish.set('Disp',Disp)

it.set_callback('forceDisp',-1)

# plot them
it.command("""
[Force = 0]
[Disp = 0]
fish history @Force
fish history @Disp
""")

it.command("model solve time 0.5")