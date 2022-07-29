# # [Getting Started](@id guide)
#
#md # ```@setup guide
#md # # this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
#md # using Plots
#md # pyplot()
#md # ```
#
#
# In this guide we introduce you to the basic functionality of this package in a step 
# by step manner.  This is a good starting point for learning about how to use this package.  
#
#-
#md # !!! tip
#md #     This guide is also available as a Jupyter notebook:
#md #     [`guide.ipynb`](@__NBVIEWER_ROOT_URL__/examples/guide.ipynb).
#-
#
# If you haven't yet, now would be a good time to install GXBeam.  It can be 
# installed from the Julia REPL by typing `]` (to enter the package manager) and then 
# running `add GXBeam`.

# Now, that the package is installed we need to load it so that we can use it.  It's also 
# often helpful to load the LinearAlgebra package.

using GXBeam, LinearAlgebra
#!jl nothing #hide

# The geometry we will be working with is a rotating beam with a swept tip as pictured.
# 
# ![](../assets/rotating-drawing.svg)
# 
# This geometry has a fixed boundary condition on the left side of the beam and rotates 
# around a point 2.5 inches to the left of the beam.  We will investigate the steady 
# behavior of this system for a variety of rotation rates when the sweep angle is 45°.
#
# ## Creating an Assembly
# 
# The first step for any analysis is to create an object of type [`Assembly`](@ref).  This 
# object stores the properties of each of the points and beam elements in our model.  
#
# To create an object of type Assembly we need the following:
#  - An array of points
#  - The starting point for each beam element
#  - The ending point for each beam element
#  - The frame of reference for each beam element, specified as a 3x3 direction cosine matrix
#  - The stiffness or compliance matrix for each beam element
#  - The mass or inverse mass matrix for each beam element, for dynamic simulations
#  - The element length and midpoint, if the element is curved
#
# We will first focus on the geometry.  We start by defining the straight section of the 
# beam.  This section extends from (2.5, 0, 0) to (34, 0, 0).  The local coordinate frame 
# for this section of the beam is the same as the global coordinate frame.  We will 
# discretize this section into 10 elements.
#
# To aid with constructing the geometry we can use the [`discretize_beam`](@ref) function.  
# We pass in the length, starting point, and number of elements of the beam section to the 
# [`discretize_beam`](@ref) function.  The function returns the lengths, endpoints, 
# midpoints, and reference frame of each beam element.

## straight section of the beam
L_b1 = 31.5 # length of straight section of the beam in inches
r_b1 = [2.5, 0, 0] # starting point of straight section of the beam
nelem_b1 = 10 # number of elements in the straight section of the beam
lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, nelem_b1)
#!jl nothing #hide

# The length of each beam element is equal since we used the number of elements to define 
# the discretization.  For more control over the discretization we can pass in the 
# discretization directly.  The following is an equally valid method for obtaining 
# uniformly spaced beam elements.

## normalized element endpoints in the straight section of the beam
disc_b1 = range(0, 1, length=nelem_b1+1)

## discretize straight beam section
lengths_b1, xp_b1, xm_b1, Cab_b1 = discretize_beam(L_b1, r_b1, disc_b1)

#!jl nothing #hide

# We now create the geometry for the swept portion of the wing.  To do so we use the same 
# [`discretize_beam`](@ref) function, but use the additional keyword argument `frame` in 
# order to define the undeformed local beam frame.  The direction cosine matrix which 
# describes the local beam frame is
# ```math
# \begin{bmatrix}
# e_{1,x} & e_{2,x} & e_{3,x} \\
# e_{1,y} & e_{2,y} & e_{3,y} \\
# e_{1,z} & e_{2,z} & e_{3,z} \\
# \end{bmatrix}
# ```
# where ``e_1``, ``e_2``, and ``e_3`` are unit vectors which define 
# the axes of the local frame of reference in the body frame of reference.  This matrix 
# may be interpreted as a transformation matrix from the undeformed local beam frame to 
# the body frame.

sweep = 45 * pi/180

## swept section of the beam
L_b2 = 6 # length of swept section of the beam
r_b2 = [34, 0, 0] # starting point of swept section of the beam
nelem_b2 = 5 # number of elements in swept section of the beam
e1 = [cos(sweep), -sin(sweep), 0] # axis 1
e2 = [sin(sweep), cos(sweep), 0] # axis 2
e3 = [0, 0, 1] # axis 3
frame_b2 = hcat(e1, e2, e3) # transformation matrix from local to body frame
lengths_b2, xp_b2, xm_b2, Cab_b2 = discretize_beam(L_b2, r_b2, nelem_b2;
    frame = frame_b2)
#!jl nothing #hide

# We will now manually combine the results of our two calls to `discretize_beam`.  Since 
# the last endpoint from the straight section is the same as the first endpoint of the 
# swept section we drop one of the endpoints when combining our results.

## combine elements and points into one array
nelem = nelem_b1 + nelem_b2 # total number of elements
points = vcat(xp_b1, xp_b2[2:end]) # all points in our assembly
start = 1:nelem_b1 + nelem_b2 # starting point of each beam element in our assembly
stop = 2:nelem_b1 + nelem_b2 + 1 # ending point of each beam element in our assembly
lengths = vcat(lengths_b1, lengths_b2) # length of each beam element in our assembly
midpoints = vcat(xm_b1, xm_b2) # midpoint of each beam element in our assembly
Cab = vcat(Cab_b1, Cab_b2) # transformation matrix from local to body frame
                           ## for each beam element in our assembly
#!jl nothing #hide

# Next we need to define the stiffness (or compliance) and mass matrices for each 
# beam element.
# 
# The compliance matrix is defined by the following equation
# ```math
# \begin{bmatrix}
# \gamma_{11} \\
# 2\gamma_{12} \\
# 2\gamma_{13} \\
# \kappa_{1} \\
# \kappa_{2} \\
# \kappa_{3}
# \end{bmatrix}
# =
# \begin{bmatrix}
#    S_{11} & S_{12} & S_{13} & S_{14} & S_{15} & S_{16} \\
#    S_{12} & S_{22} & S_{23} & S_{24} & S_{25} & S_{26} \\
#    S_{13} & S_{23} & S_{33} & S_{34} & S_{35} & S_{36} \\
#    S_{14} & S_{24} & S_{43} & S_{44} & S_{45} & S_{46} \\
#    S_{15} & S_{25} & S_{35} & S_{45} & S_{55} & S_{56} \\
#    S_{16} & S_{26} & S_{36} & S_{46} & S_{56} & S_{66}
# \end{bmatrix}
# \begin{bmatrix}
#    F_{1} \\
#    F_{2} \\
#    F_{3} \\
#    M_{1} \\
#    M_{2} \\
#    M_{3}
# \end{bmatrix}
# ```
# with the variables defined as follows:
#  - ``\gamma_{11}``: beam axial strain
#  - ``2\gamma_{12}`` engineering transverse strain along axis 2
#  - ``2\gamma_{13}`` engineering transverse strain along axis 3
#  - ``\kappa_1``: twist
#  - ``\kappa_2``: curvature about axis 2
#  - ``\kappa_3``: curvature about axis 3
#  - ``F_i``: resultant force about axis i
#  - ``M_i``: resultant moment about axis i
#
# The mass matrix is defined using the following equation
# ```math
# \begin{bmatrix}
#    P_{1} \\
#    P_{2} \\
#    P_{3} \\
#    H_{1} \\
#    H_{2} \\
#    H_{3}
# \end{bmatrix}
# =
# \begin{bmatrix}
#    \mu & 0 & 0 & 0 & \mu x_{m3} & -\mu x_{m2} \\
#    0 & \mu & 0 & -\mu x_{m3} & 0 & 0 \\
#    0 & 0 & \mu & \mu x_{m2} & 0 & 0 \\
#    0 & -\mu x_{m3} & \mu x_{m2} & i_{22} + i_{33} & 0 & 0 \\
#    \mu x_{m3}  & 0 & 0 & 0 & i_{22} & -i_{23} \\
#    -\mu x_{m2} & 0 & 0 & 0 & -i_{23} & i_{33}
# \end{bmatrix}
# \begin{bmatrix}
#    V_{1} \\
#    V_{2} \\
#    V_{3} \\
#    \Omega_{1} \\
#    \Omega_{2} \\
#    \Omega_{3}
# \end{bmatrix}
# ```
# with the variables defined as follows:
#  - ``P``: linear momentum per unit length
#  - ``H``: angular momentum per unit length
#  - ``V``: linear velocity
#  - ``\Omega``: angular velocity
#  - ``\mu``: mass per unit length
#  - ``(x_{m2}, x_{m3})``: mass center location
#  - ``i_{22}``: mass moment of inertia about axis 2
#  - ``i_{33}``: mass moment of inertia about axis 3
#  - ``i_{23}``: product of inertia
#
# We assume that our beam has a constant cross section with the following properties:
#  - 1 inch width
#  - 0.063 inch height
#  - 1.06 x 10^7 lb/in^2 elastic modulus
#  - 0.325 Poisson's ratio
#  - 2.51 x 10^-4 lb sec^2/in^4 density
#
# We also assume the following shear and torsion correction factors:
#  - ``k_y = 1.2000001839588001``
#  - ``k_z = 14.625127919304001``
#  - ``k_t = 65.85255016982444``

## cross section
w = 1 # inch
h = 0.063 # inch

## material properties
E = 1.06e7 # lb/in^2
ν = 0.325
ρ = 2.51e-4 # lb sec^2/in^4

## shear and torsion correction factors
ky = 1.2000001839588001
kz = 14.625127919304001
kt = 65.85255016982444

A = h*w
Iyy = w*h^3/12
Izz = w^3*h/12
J = Iyy + Izz

## apply corrections
Ay = A/ky
Az = A/kz
Jx = J/kt

G = E/(2*(1+ν))

compliance = fill(Diagonal([1/(E*A), 1/(G*Ay), 1/(G*Az), 1/(G*Jx), 1/(E*Iyy),
    1/(E*Izz)]), nelem)

mass = fill(Diagonal([ρ*A, ρ*A, ρ*A, ρ*J, ρ*Iyy, ρ*Izz]), nelem)
#!jl nothing #hide

# Our case is simple enough that we can analytically calculate most values for the 
# compliance and mass matrices, but this is not generally the case.  For more complex 
# geometries/structures see the section below titled [`Section Properties`](@ref)
# Also note that any row/column of the stiffness and/or compliance matrix which is zero 
# will be interpreted as infinitely stiff in that degree of freedom.  This corresponds to a 
# row/column of zeros in the compliance matrix.
#
# We are now ready to put together our assembly.

assembly = Assembly(points, start, stop;
   compliance = compliance,
   mass = mass,
   frames = Cab,
   lengths = lengths,
   midpoints = midpoints)
#!jl nothing #hide

# At this point this is probably a good time to check that the geometry of our assembly 
# is correct.  We can do this by visualizing the geometry in 
# [ParaView](https://www.paraview.org/).  We can use the [`write_vtk`](@ref) function to 
# do this.  Note that in order to visualize the generated file yourself you will need to 
# [install ParaView](https://www.paraview.org/download/) separately.

write_vtk("rotating-geometry", assembly)

# ![](../assets/rotating-geometry.png)

# ## Section Properties

# In the above example, and in most of the examples that follow, the stiffness and inertial properties 
# are provided as direct inputs. In general, however, we need to compute these section properties from 
# input geometry and material properties.  

# We implemented the 2D finite element procedure described by 
# [Giavotto et al.](https://doi.org/10.1016/0045-7949(83)90179-7) and that is also described 
# in the [BECAS User Guide](https://backend.orbit.dtu.dk/ws/portalfiles/portal/7711204/ris_r_1785.pdf).  
# This approach is much more accurate than those based on classical laminate theory, especially for coupled 
# and transverse stiffnesses, and of course is much more accurate for thicker sections.  
# Our implementation uses bilinear quadrilateral (Q4) elements, and is tailored 
# for use in design optimization.  First, the implementation is fast.  The memory for all large matrices is 
# preallocated, in sparse formats where applicable.  This allows operations to be done in place for subsequent 
# iterations within an optimization. We also save the matrix factorization in the linear solve and use sparse 
# solvers.   Second, meshes are carefully resized to avoid discrete jumps in properties as airfoil/thickness 
# dimensions change.

# There are two main ways to compute section properties using this package.  The first is from explicit creation 
# of a mesh, i.e., nodes and elements. This is the most general approach, and can be used for any shape.  
# However, as our primary usage is for airfoils, we have a second convenience function that generates the nodes 
# and elements for a parameterized airfoil.  This approach does not require understanding the format of nodes 
# and elements, but rather the specificiations for the paramterized airfoil.  Other convenience functions could be 
# created for other common geometries.  These two different approaches are explained below, but before describing 
# them we note the material definition as that is needed in both methods.

# #### Material Properties

# We can specify a general orthotropic [`Material`](@ref) with three independent Young's moduli, shear moduli, Poisson's 
# ratios, and the material density.  ``E_i`` is the Young's modulus along axis ``i``, ``G_{ij}`` is the shear 
# modulus in direction ``j`` on the plane whose normal is in direction ``i``, and ``\nu_{ij}`` is the Poisson’s 
# ratio that corresponds to a contraction in direction ``j`` when an extension is applied in direction ``i`` 
# (from symmetry the opposite definitions would also apply, where we swap ``i`` and ``j``).  
# Also from symmetry we know that ``\nu_{ij} E_j = \nu_{ji} E_i``.  For example:

E1 = 10e9
E2 = 5e9
E3 = 5e9
G12 = 8e9
G13 = 4e9
G23 = 4e9
nu12 = 0.3
nu13 = 0.3
nu23 = 0.3
rho = 1.8e3

mat = Material(E1, E2, E3, G12, G13, G23, nu12, nu13, nu23, rho)
#!jl nothing #hide

# The ply, element, and beam coordinate systems are shown below.  Axis 1 is along the ply main fiber direction, 
# 2 is the transverse direction, and 3 is normal to the ply. For a wing/blade 1 is nominally along the span/radius, 
# 2 is tangent to airfoil, and 3 is normal to airfoil (these are nominal as they would only be exact if the fiber 
# orientation, ``\theta``, is 0 degrees).

# ![](../assets/material-coordinate.png)

# #### Nodes and Elements

# The structure is discretized into nodes.  A [`Node`](@ref) is a point given by an x, y coordinate pair.  

x = 0.0
y = 0.0
node = Node(x, y)
#!jl nothing #hide

# Of course, a structure has many nodes, and these should be assembled in a vector of nodes.  The node number  
# is given by the vector index.  In this example we create four nodes in the shape of a square from which we can 
# later construct an element.

x = 0.0
y = 0.0
node1 = Node(x, y)
x = 1.0
y = 0.0
node2 = Node(x, y)
x = 1.0
y = 1.0
node3 = Node(x, y)
x = 0.0
y = 1.0
node4 = Node(x, y)

nodes = [node1; node2; node3; node4]
#!jl nothing #hide

# The structure is constructed with elements.  Each [`MeshElement`](@ref) is made of four nodes, ordered as shown below Each 
# element also has a material, and that material has a fiber orientation ``\theta`` (see previous figure).

# ![](../assets/element.png)

nodenumbers = [1, 2, 3, 4]
material = mat
theta = 20*pi/180
element = MeshElement(nodenumbers, material, theta)
#!jl nothing #hide

# A mesh would then consist of an array of nodes and an array of elements. As an example we construct a mesh 
# for the square cross section with isotropic material shown in the BECAS User Guide (which is just a square of 
# side length 0.1).  We discretize the square into 100 equally-sized elements, 10 in each dimension.

iso = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1000.0)
x = range(-0.05, 0.05, length=11)
y = range(-0.05, 0.05, length=11)

nodes = Vector{Node}(undef, 11*11)
elements = Vector{MeshElement}(undef, 10*10)

let 
    m = 1
    for i = 1:11
        for j = 1:11
            nodes[m] = Node(x[i], y[j])
            m += 1
        end
    end

    m = 1
    for i = 1:10
        for j = 1:10
            elements[m] = MeshElement([11*(i-1)+j, 11*(i)+j, 11*(i)+j+1, 11*(i-1)+j+1], iso, 0.0)
            m += 1
        end
    end
end
#!jl nothing #hide

# With the mesh constructed we can now compute the compliance matrix.  The [`compliance_matrix`](@ref) is computed about 
# the shear center.

S, sc, tc = compliance_matrix(nodes, elements)
#!jl nothing #hide

# S is the compliance matrix, sc is the x and y coordinates for the shear center, and tc is the x and y 
# coordinates for the tension center.  These coordinates are relative to the origin (0, 0) of the mesh.

# Inertial properties can be computed from the [`mass_matrix`](@ref) function:

M, mc = mass_matrix(nodes, elements)
#!jl nothing #hide

# where M is the mass matrix and mc is the x and y coordinates for the mass center.

# #### Airfoil Shape

# The airfoil outer mold line is defined by a set of x, y coordinates, normalized by chord. The points must 
# start at the trailing edge, traverse counterclockwise (i.e., upper surface first), and end at the trailing 
# edge as shown below.  
# The trailing edge can be blunt or sharp. For the former the trailing edge would start and end at different points, 
# and in the later they would be the same point.

# ![](../assets/airfoil.png)

# In this example we will use ST1 for the MH-104 airfoil described in the paper by 
# [Chen, Yu, and Capellaro](https://doi.org/10.1002/we.372).

xaf = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
yaf = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]
#!jl nothing #hide

# Next, we define the chord length, the twist angle, and the pitch axis.  The pitch axis is the x-coordinate, 
# measured back from the leading edge, normalized by chord, that defines where the airfoil should be twisted about
# as seen in the figure below.  
# So a pitch axis of 0.25 means that the airfoil rotates about the quarter chord.  Positive twist is in the 
# direction of increasing angle of attack. In our case the twist is zero so the pitch axis is irrelevant though 
# we use the value noted in the paper.

# ![](../assets/twist.png)

chord = 1.9
twist = 0.0*pi/180
paxis = 0.4750 / chord
#!jl nothing #hide

# We now define the materials we will use in our layup.  There are five materials in this case, and for 
# convenience we put them in an array so we can refer to them by number.

uni = Material(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
double = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
gelcoat = Material(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
nexus = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
balsa = Material(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
mat = [uni, double, gelcoat, nexus, balsa]
#!jl nothing #hide

# The heart of the parmeterization is defining the ply sequence.  A single ply (or multiple plys of the same 
# material and orientation) is given by a Layer.  A Layer is defined by a material, a thickness, and a fiber 
# orientation. Thicknesses are absolute (not normalized), and orientation angles should be 
# in radians.  For example:

t = 0.001
theta = 20*pi/180
layer = Layer(balsa, t, theta)
#!jl nothing #hide

# A given segment on the airfoil will have multiple layers, which defines the ply sequencing at that section
# (see figure below).  
# ![](../assets/layup.png)
# For example the first segment of this airfoil uses three materials: a gelcoat followed by nexus then the 
# double-bias FRP.  Note that the sequencing always starts from the outer edge of the airfoil.  Recall that we 
# placed our materials in an array so we can refer to them by number.  We then use broadcasting to create a 
# vector of these three layers.

idx = [3, 4, 2]  # material index
t = [0.000381, 0.00051, 18*0.00053]
theta = [0, 0, 20]*pi/180
layup1 = Layer.(mat[idx], t, theta)
#!jl nothing #hide

# This particular airfoil is made up of four separate segments, so we now define the next three using the 
# specifications from the document.

idx = [3, 4, 2]
t = [0.000381, 0.00051, 33*0.00053]
theta = [0, 0, 20]*pi/180
layup2 = Layer.(mat[idx], t, theta)

idx = [3, 4, 2, 1, 5, 1, 2]
t = [0.000381, 0.00051, 17*0.00053, 38*0.00053, 1*0.003125, 37*0.00053, 16*0.00053]
theta = [0, 0, 20, 30, 0, 30, 20]*pi/180
layup3 = Layer.(mat[idx], t, theta)

idx = [3, 4, 2, 5, 2]
t = [0.000381, 0.00051, 17*0.00053, 0.003125, 16*0.00053]
theta = [0, 0, 20, 0, 0]*pi/180
layup4 = Layer.(mat[idx], t, theta) 
#!jl nothing #hide

# These four layups correspond four different regions of the airfoil.  We concatenate them into a vector 
# called segments.

segments = [layup1, layup2, layup3, layup4]
#!jl nothing #hide

# We need to define over what region each layup corresponds to.  We do this by providing the normalized x 
# locations as breakpoints.  

xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]
#!jl nothing #hide

# The above means that layup 1 applies from x/c = 0 to 0.0041, layup 2 from x/c = 0.0041 to 0.1147, and so on.

# Shear webs are defined in the same way as segments. Each web has a stack of layers, and the airfoil can have 
# as many webs as desired.  In our case, we have two webs, both with the same ply stack.  The ordering if from 
# the leading edge side towards the trailing edge, although usually the web layups are symmetric anyway.

idx = [1, 5, 1]
t = [38*0.00053, 0.003125, 38*0.00053]
theta = [0, 0, 0]*pi/180
web = Layer.(mat[idx], t, theta)

webs = [web, web]
#!jl nothing #hide

# We also need to define the normalized x location of each web center.
webloc = [0.15, 0.5]
#!jl nothing #hide
# Thus, our first web is at x/c = 0.15, the second at x/c = 0.5.

# We now have all the information we need to create the structural mesh using the [`afmesh`](@ref) function.
nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)
#!jl nothing #hide

# We can now compute [`compliance_matrix`](@ref) and inertial matices ([`mass_matrix`](@ref)) in the same 
# manner as shown in the prior section.
S, sc, tc = compliance_matrix(nodes, elements)
M, mc = mass_matrix(nodes, elements)
#!jl nothing #hide

# For optimization usage, one should also initialize the cache upfront (see [`initialize_cache`](@ref)).  
# Then re-use that cache for subsequent 
# calls to `compliance_matrix`.  That wiill avoid recreating the cache at each iteration and will speed up computation.  
# The cache does not needed to be re-initialized as long as the matrix sizes stay constant 
# (the geometry can change sizes, but the overall mesh connectivity remains the same, which should be the case
# during gradient-based design optimization)
cache = initialize_cache(nodes, elements)
S, sc, tc = compliance_matrix(nodes, elements, cache=cache)  # these calls are now much faster
S, sc, tc = compliance_matrix(nodes, elements, cache=cache)
#!jl nothing #hide

# For visualizing the mesh there is a provided helper function [`plotmesh`](@ref).
# It uses PyPlot, although that must be loaded on 
# the user-side and passed in as it is not a dependency in the package.

# ```julia
# using PyPlot
# figure(); plotmesh(node, elements, PyPlot)
# ```

# ![](../assets/afmesh.png)



# #### Airfoil Mesh Control

# A few additional parameters are available to control the airfoil meshing.  By default the mesh is discretized 
# in the tangential direction using the airfoil coordinates. In other words, to get a finer mesh in this direction 
# we could specify more airfoil points.  By default the mesh in the normal direction uses the layers defined 
# by the segments.  Since we use the same number of elements for each segment, we subdivide some layers as needed 
# (and thus the segment with the most layers defines the number of elements).  In our above example the the 
# third segment has the most layers (7) so we discretize the mesh with 7 elements in the tangential direction 
# for all segments.  To get a finer mesh in this direction we could breakup a segment into more layers 
# (for example a layer of material 1 of 0.1 thickness could be divided into four layers of material 1 with 
# 0.025 thickness).

# However, doing this remeshing manually is a bit tedious so there are functions that can automate this.  
# First, there are two optional parameters that control the desired spacing in the tangential and normal 
# directions.  The first is "ds", which if provided, will resample the airfoil with approximately this spacing, 
# normalized by chord.  For example, ds=0.01 will create points on the airfoil spaced roughly 1% chord apart.  
# Siimilarly, there is an optional parameter "dt", which resamples the thicknesses with this maximium mesh size 
# (thickness is absolute). For example, dt=0.01, will target a maximum element thickness of 0.01. Recall that 
# the total number of elements remains constant along the airfoil, so most thicknesses will less than this
# value. This same parameter is used for discretizing along the thickness of the webs.  

# These two parameters are simple to use, but are typically not desirable for gradient-based optimization if 
# ply thicknesses or airfoil coordinates are changed.  This is because a fixed element size will experience 
# discrete jumps as the airfoil is resized.  Thus, rather than specify a fixed element size, we should specify 
# a fixed number of elements.  The optional parameters ns and nt provide this functionality.  

# The parameter ns is an array specifying the number of elements to create in each segment.  In our case we have 
# four different segments so we could define an array like: `ns = [15, 20, 40, 30]`.  This is less convenient 
# as it is not as immediately obvious how many elements should to put in each sector in order to provide a 
# consistently spaced mesh.  However, as the airfoil shape changes, the mesh will simply stretch/shrink and 
# never experience discrete jumps.  Similarly, we can specify nt to define the number of elements to use 
# across the thickness.  The parameter nt is actually an array of arrays.  Each subaray defines how many 
# elements should be used for each layer.  For example 
# `nt = [[1, 1, 5], [1, 1, 5], [1, 1, 1, 1, 1, 1, 1], [1, 1, 2, 1, 2]]` would subdivide the innermost 
# layer of the first segment into 5 elements, and so on.

# Finally, for discretizing along the webs we can use `wns` and `wnt`.  The parmeter wns is an integer and 
# specifies the number of elements to use along the height of the web (same definition as ns but a web only 
# needs one number). The default is four. The parameters wnt has the same definition as nt but applies to the webs.

# ## Point Masses
#
# We won't be applying point masses to our model, but we will demonstrate how to do so.
# 
# Point masses are defined by using the constructor [`PointMass`](@ref) and may be attached
# to any point.  One instance of [`PointMass`](@ref) must be created for every point 
# with attached point masses.  These instances of [`PointMass`](@ref) are then stored 
# in a dictionary with keys corresponding to each point index.
# 
# Each [`PointMass`](@ref) contains a 6x6 mass matrix which describes the relationship 
# between the linear/angular velocity of the point and the linear/angular momentum 
# of the point mass.  For a single point mass, this matrix is defined as
# ```math
# \begin{bmatrix}
#    P_{x} \\
#    P_{y} \\
#    P_{z} \\
#    H_{x} \\
#    H_{y} \\
#    H_{z}
# \end{bmatrix}
# =
# \begin{bmatrix}
#    m & 0 & 0 & 0 & m p_{z} & -m p_{y} \\
#    0 & m & 0 & -m p_{z} & 0 & m p_{x} \\
#    0 & 0 & m & m p_{y} & -m p_{x} & 0 \\
#    0 & -m p_{z} & m p_{y} & I_{xx}^* & -I_{xy}^* & -I_{xz}^* \\
#    m p_{z}  & 0 & -m p_{x} & -I_{xy}^* & I_{yy}^* & -I_{yz}^* \\
#    -m p_{y} & m p_{x} & 0 & -I_{xz}^* & -I_{yz}^* & I_{zz}^*
# \end{bmatrix}
# \begin{bmatrix}
#    V_{x} \\
#    V_{y} \\
#    V_{z} \\
#    \Omega_{x} \\
#    \Omega_{y} \\
#    \Omega_{z}
# \end{bmatrix}
# ```
# where ``m`` is the mass of the point mass, ``p`` is the position of the point mass 
# relative to the point to which it is attached, and ``I^*`` is the 
# inertia matrix corresponding to the point mass, defined relative to the point.  
# Multiple point masses may be modeled by adding their respective mass 
# matrices together.

# Objects of type [`PointMass`](@ref) may be constructed by providing the fully populated 
# mass matrix as described above or by providing the mass, offset, and inertia matrix of 
# the point mass, with the later being the inertia matrix of the point mass about its 
# center of gravity rather than the beam center.  To demonstrate, the following code places 
# a 10 kg tip mass at the end of our swept beam.

m = 10 # mass
p = zeros(3) # relative location
J = zeros(3,3) # inertia matrix (about the point mass center of gravity)

## create dictionary of point masses
point_masses = Dict(
    nelem+1 => PointMass(m, p, J)
    )

#!jl nothing #hide

# ## Defining Distributed Loads
# 
# We won't be applying distributed loads to our model, but we will demonstrate how to do so.
# 
# Distributed loads are defined by using the constructor [`DistributedLoads`](@ref).  One 
# instance of [`DistributedLoads`](@ref) must be created for every beam element on which 
# the distributed load is applied.  These instances of [`DistributedLoads`](@ref) are then 
# stored in a dictionary with keys corresponding to each beam element index.  
#
# To define a [`DistributedLoads`](@ref) the assembly, element number, and distributed 
# load functions must be passed to [`DistributedLoads`](@ref).  Possible distributed 
# load functions are:
# - `fx`: Distributed x-direction force
# - `fy`: Distributed y-direction force
# - `fz`: Distributed z-direction force
# - `mx`: Distributed x-direction moment
# - `my`: Distributed y-direction moment
# - `mz`: Distributed z-direction moment
# - `fx_follower`: Distributed x-direction follower force
# - `fy_follower`: Distributed y-direction follower force
# - `fz_follower`: Distributed z-direction follower force
# - `mx_follower`: Distributed x-direction follower moment
# - `my_follower`: Distributed y-direction follower moment
# - `mz_follower`: Distributed z-direction follower moment
# 
# Each of these forces/moments are specified as functions of the arbitrary coordinate ``s```.  
# The ``s``-coordinate at the start and end of the beam element can be specified 
# using the keyword arguments ``s1`` and ``s2``.

# For example, the following code applies a uniform 10 pound distributed load in the 
# z-direction on all beam elements:

distributed_loads = Dict{Int64, DistributedLoads{Float64}}()
for ielem in 1:nelem
    distributed_loads[ielem] = DistributedLoads(assembly, ielem; fz = (s) -> 10)
end
#!jl nothing #hide

# To instead use a follower force (a force that rotates with the structure) we would use 
# the following code:

distributed_loads = Dict{Int64, DistributedLoads{Float64}}()
for ielem in 1:nelem
    distributed_loads[ielem] = DistributedLoads(assembly, ielem;
        fz_follower = (s) -> 10)
end
#!jl nothing #hide

# The units are arbitrary, but must be consistent with the units used when constructing 
# the beam assembly.  Also note that both non-follower and follower forces may exist 
# simultaneously.
# 
# Note that the distributed loads are integrated over each element when they 
# are created using 4-point Gauss-Legendre quadrature.  If more control over the 
# integration is desired one may specify a custom integration method as described in the 
# documentation for [`DistributedLoads`](@ref).

# ## Defining Prescribed Conditions
# 
# Whereas distributed loads are applied to beam elements, prescribed conditions are 
# external loads or displacement boundary conditions applied to points. One instance of 
# [`PrescribedConditions`](@ref) must be created for every point on which prescribed 
# conditions are applied.  These instances of `PrescribedConditions` are then stored in a 
# dictionary with keys corresponding to each point index.
# 
# Possible prescribed conditions include:
# - `ux`: Prescribed x-displacement
# - `uy`: Prescribed y-displacement
# - `uz`: Prescribed z-displacement
# - `theta_x`: Prescribed first Wiener-Milenkovic parameter
# - `theta_y`: Prescribed second Wiener-Milenkovic parameter
# - `theta_z`: Prescribed third Wiener-Milenkovic parameter
# - `Fx`: Prescribed x-direction force
# - `Fy`: Prescribed y-direction force
# - `Fz`: Prescribed z-direction force
# - `Mx`: Prescribed x-axis moment
# - `My`: Prescribed y-axis moment
# - `Mz`: Prescribed z-axis moment
# - `Fx_follower`: Prescribed x-direction follower force
# - `Fy_follower`: Prescribed y-direction follower force
# - `Fz_follower`: Prescribed z-direction follower force
# - `Mx_follower`: Prescribed x-direction follower moment
# - `My_follower`: Prescribed y-direction follower moment
# - `Mz_follower`: Prescribed z-direction follower moment
# 
# One can apply both force and displacement boundary conditions to the same point, but one 
# cannot specify a force and displacement condition at the same point corresponding 
# to the same degree of freedom.
# 
# Here we create a fixed boundary condition on the left side of the beam.

## create dictionary of prescribed conditions
prescribed_conditions = Dict(
    ## root section is fixed
    1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

#!jl nothing #hide

# Note that most problems should have at least one point where deflections and/or 
# rotations are constrained in order to be well-posed.

# ## Pre-Allocating Memory for an Analysis
# 
# At this point we have everything we need to perform an analysis.  However, since we will 
# be performing multiple analyses using the same assembly we can save computational time 
# be pre-allocating memory for the analysis.  This can be done by constructing an object of
# type [`AbstractSystem`](@ref).  There are two main options: [`StaticSystem`](@ref) for 
# static systems and [`DynamicSystem`](@ref) for dynamic systems.  The third option:
# [`ExpandedSystem`](@ref) is primarily useful when constructing a constant mass matrix 
# system for use with [`DifferentialEquations`](https://github.com/SciML/DifferentialEquations.jl)  
# Since our system is rotating, we construct an object of type [`DynamicSystem`](@ref).

system = DynamicSystem(assembly)
#!jl nothing #hide

# ## Performing a Steady State Analysis
# 
# We're now ready to perform our steady state analyses.  This can be done by calling 
# [`steady_state_analysis!`](@ref) with the pre-allocated system storage, assembly, 
# angular velocity, and the prescribed point conditions.  A linear analysis may be 
# performed instead of a nonlinear analysis by using the `linear` keyword argument.
# 
# After each analysis we'll also construct an object of type [`AssemblyState`](@ref) so 
# that we can save the results of each analysis prior to re-using the pre-allocated 
# memory for the next analysis.

rpm = 0:25:750

linear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
for i = 1:length(rpm)

    ## global frame rotation
    w0 = [0, 0, rpm[i]*(2*pi)/60]

    ## perform linear steady state analysis
    _, converged = steady_state_analysis!(system, assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions,
        linear = true)

    linear_states[i] = AssemblyState(system, assembly;
        prescribed_conditions=prescribed_conditions)

end

reset_state!(system)

nonlinear_states = Vector{AssemblyState{Float64}}(undef, length(rpm))
for i = 1:length(rpm)

   ## global frame rotation
   w0 = [0, 0, rpm[i]*(2*pi)/60]

    ## perform nonlinear steady state analysis
    _, converged = steady_state_analysis!(system, assembly,
        angular_velocity = w0,
        prescribed_conditions = prescribed_conditions)

     nonlinear_states[i] = AssemblyState(system, assembly;
         prescribed_conditions=prescribed_conditions)

end

#!jl nothing #hide

# ## Post Processing Results
#
# We can access the fields in each instance of [`AssemblyState`](@ref) in order to plot 
# various quantities of interest.  This object stores an array of objects of type 
# [`PointState`](@ref) in the field `points` and an array of objects of type 
# [`ElementState`](@ref) in the field `elements`.  
# 
# The fields of [`PointState`](@ref) are the following:
#  - `u`: point linear displacement (in the global frame)
#  - `theta`: point angular displacement (in the global frame)
#  - `F`: externally applied forces on the point (in the global frame)
#  - `M`: externally applied moments on the point (in the global frame)
#  - `V`: linear velocity (in the global frame)
#  - `Omega`: angular velocity (in the global frame)
# 
# The fields of [`ElementState`](@ref) are the following:
#  - `u`: element displacement (in the global frame )
#  - `theta`: angular displacement (in the global frame)
#  - `F`: resultant forces (in the local element frame)
#  - `M`: resultant moments (in the local element frame)
#  - `V`: linear velocity (in the global frame)
#  - `Omega`: angular velocity (in the global frame)
# 
# Angular displacements are expressed in terms of Wiener-Milenkovic parameters. 
#
# To demonstrate how these fields can be accessed we will now plot the root moment and 
# tip deflections.

using Plots
#md using Suppressor #hide
pyplot()
#!jl nothing #hide

#-

#md @suppress_err begin #hide

## root moment
plot(
    xlim = (0, 760),
    xticks = 0:100:750,
    xlabel = "Angular Speed (RPM)",
    yticks = 0.0:2:12,
    ylabel = "\$M_z\$ at the root (lb-in)",
    grid = false,
    overwrite_figure=false
    )
Mz_nl = [-nonlinear_states[i].points[1].M[3] for i = 1:length(rpm)]
Mz_l = [-linear_states[i].points[1].M[3] for i = 1:length(rpm)]
plot!(rpm, Mz_nl, label="Nonlinear")
plot!(rpm, Mz_l, label="Linear")
plot!(show=true) #!nb
#md savefig("../assets/guide-Mz.svg") #hide
#md closeall() #hide
#md end #hide
#md nothing #hide

#md # ![](../assets/guide-Mz.svg)

#-

#md @suppress_err begin #hide

## x tip deflection
plot(
    xlim = (0, 760),
    xticks = 0:100:750,
    xlabel = "Angular Speed (RPM)",
    ylim = (-0.002, 0.074),
    yticks = 0.0:0.01:0.07,
    ylabel = "\$u_x\$ at the tip (in)",
    grid = false,
    overwrite_figure=false
    )
ux_nl = [nonlinear_states[i].points[end].u[1] for i = 1:length(rpm)]
ux_l = [linear_states[i].points[end].u[1] for i = 1:length(rpm)]
plot!(rpm, ux_nl, label="Nonlinear")
plot!(rpm, ux_l, label="Linear")
plot!(show=true) #!nb
#md savefig("../assets/guide-ux.svg") #hide
#md closeall() #hide
#md end #hide
#md nothing #hide

#md # ![](../assets/guide-ux.svg)

#-

#md @suppress_err begin #hide

## y tip deflection
plot(
    xlim = (0, 760),
    xticks = 0:100:750,
    xlabel = "Angular Speed (RPM)",
    ylim = (-0.01, 0.27),
    yticks = 0.0:0.05:0.25,
    ylabel = "\$u_y\$ at the tip (in)",
    grid = false,
    overwrite_figure=false
    )
uy_nl = [nonlinear_states[i].points[end].u[2] for i = 1:length(rpm)]
uy_l = [linear_states[i].points[end].u[2] for i = 1:length(rpm)]
plot!(rpm, uy_nl, label="Nonlinear")
plot!(rpm, uy_l, label="Linear")
plot!(show=true) #!nb
#md savefig("../assets/guide-uy.svg") #hide
#md closeall() #hide
#md end #hide
#md nothing #hide

#md # ![](../assets/guide-uy.svg)

#-

#md @suppress_err begin #hide

## rotation of the tip
plot(
    xlim = (0, 760),
    xticks = 0:100:750,
    xlabel = "Angular Speed (RPM)",
    ylabel = "\$θ_z\$ at the tip",
    grid = false,
    overwrite_figure=false
    )
theta_z_nl = [4*atan(nonlinear_states[i].points[end].theta[3]/4)
    for i = 1:length(rpm)]
theta_z_l = [4*atan(linear_states[i].points[end].theta[3]/4)
    for i = 1:length(rpm)]

plot!(rpm, theta_z_nl, label="Nonlinear")
plot!(rpm, theta_z_l, label="Linear")
plot!(show=true) #!nb
#md savefig("../assets/guide-theta_z.svg") #hide
#md closeall() #hide
#md end #hide
#md nothing #hide

#md # ![](../assets/guide-theta_z.svg)

#-

# ## Other Capabilities
# 
# Further information about how to use this package may be obtained by looking through the
# examples or browsing the [Public API](@ref).
