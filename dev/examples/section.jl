using GXBeam

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

x = 0.0
y = 0.0
node = Node(x, y)

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

nodenumbers = [1, 2, 3, 4]
material = mat
theta = 20*pi/180
element = MeshElement(nodenumbers, material, theta)

iso = Material(100.0, 100.0, 100.0, 41.667, 41.667, 41.667, 0.2, 0.2, 0.2, 1000.0)
x = range(-0.05, 0.05, length=11)
y = range(-0.05, 0.05, length=11)

nodes = Vector{Node{Float64}}(undef, 11*11)
elements = Vector{MeshElement{Vector{Int64},Float64}}(undef, 10*10)

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

S, sc, tc = compliance_matrix(nodes, elements)

M, mc = mass_matrix(nodes, elements)

xaf = [1.00000000, 0.99619582, 0.98515158, 0.96764209, 0.94421447, 0.91510964, 0.88074158, 0.84177999, 0.79894110, 0.75297076, 0.70461763, 0.65461515, 0.60366461, 0.55242353, 0.50149950, 0.45144530, 0.40276150, 0.35589801, 0.31131449, 0.26917194, 0.22927064, 0.19167283, 0.15672257, 0.12469599, 0.09585870, 0.07046974, 0.04874337, 0.03081405, 0.01681379, 0.00687971, 0.00143518, 0.00053606, 0.00006572, 0.00001249, 0.00023032, 0.00079945, 0.00170287, 0.00354717, 0.00592084, 0.01810144, 0.03471169, 0.05589286, 0.08132751, 0.11073805, 0.14391397, 0.18067874, 0.22089879, 0.26433734, 0.31062190, 0.35933893, 0.40999990, 0.46204424, 0.51483073, 0.56767889, 0.61998250, 0.67114514, 0.72054815, 0.76758733, 0.81168064, 0.85227225, 0.88883823, 0.92088961, 0.94797259, 0.96977487, 0.98607009, 0.99640466, 1.00000000]
yaf = [0.00000000, 0.00017047, 0.00100213, 0.00285474, 0.00556001, 0.00906779, 0.01357364, 0.01916802, 0.02580144, 0.03334313, 0.04158593, 0.05026338, 0.05906756, 0.06766426, 0.07571157, 0.08287416, 0.08882939, 0.09329359, 0.09592864, 0.09626763, 0.09424396, 0.09023579, 0.08451656, 0.07727756, 0.06875796, 0.05918984, 0.04880096, 0.03786904, 0.02676332, 0.01592385, 0.00647946, 0.00370956, 0.00112514, -0.00046881, -0.00191488, -0.00329201, -0.00470585, -0.00688469, -0.00912202, -0.01720842, -0.02488211, -0.03226730, -0.03908459, -0.04503763, -0.04986836, -0.05338180, -0.05551392, -0.05636585, -0.05605816, -0.05472399, -0.05254383, -0.04969990, -0.04637175, -0.04264894, -0.03859653, -0.03433153, -0.02996944, -0.02560890, -0.02134397, -0.01726049, -0.01343567, -0.00993849, -0.00679919, -0.00402321, -0.00180118, -0.00044469, 0.00000000]

chord = 1.9
twist = 0.0*pi/180
paxis = 0.4750 / chord

uni = Material(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3)
double = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3)
gelcoat = Material(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3)
nexus = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3)
balsa = Material(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3)
mat = [uni, double, gelcoat, nexus, balsa]

t = 0.001
theta = 20*pi/180
layer = Layer(balsa, t, theta)

idx = [3, 4, 2]  # material index
t = [0.000381, 0.00051, 18*0.00053]
theta = [0, 0, 20]*pi/180
layup1 = Layer.(mat[idx], t, theta)

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

segments = [layup1, layup2, layup3, layup4]

xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]

idx = [1, 5, 1]
t = [38*0.00053, 0.003125, 38*0.00053]
theta = [0, 0, 0]*pi/180
web = Layer.(mat[idx], t, theta)

webs = [web, web]

webloc = [0.15, 0.5]

nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)

S, sc, tc = compliance_matrix(nodes, elements)
M, mc = mass_matrix(nodes, elements)

cache = initialize_cache(nodes, elements)
S, sc, tc = compliance_matrix(nodes, elements, cache=cache)  # these calls are now much faster
S, sc, tc = compliance_matrix(nodes, elements, cache=cache)

F = [0.0; 0; 0]
M = [2e3; 1e3; 0]
strain_b, stress_b, strain_p, stress_p = strain_recovery(F, M, nodes, elements, cache)

uni = Material(37.00e9, 9.00e9, 9.00e9, 4.00e9, 4.00e9, 4.00e9, 0.28, 0.28, 0.28, 1.86e3,
    37.00e9/100, 37.00e9/100, 9.00e9/100, 9.00e9/100, 9.00e9/100, 9.00e9/100, 4.00e9/100, 4.00e9/100, 4.00e9/100)
double = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.83e3,
    10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 8.00e9/100, 8.00e9/100, 8.00e9/100)
gelcoat = Material(1e1, 1e1, 1e1, 1.0, 1.0, 1.0, 0.30, 0.30, 0.30, 1.83e3,
    1e1/100, 1e1/100, 1e1/100, 1e1/100, 1e1/100, 1e1/100, 1.0/100, 1.0/100, 1.0/100)
nexus = Material(10.30e9, 10.30e9, 10.30e9, 8.00e9, 8.00e9, 8.00e9, 0.30, 0.30, 0.30, 1.664e3,
    10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 10.30e9/100, 8.00e9/100, 8.00e9/100, 8.00e9/100)
balsa = Material(0.01e9, 0.01e9, 0.01e9, 2e5, 2e5, 2e5, 0.30, 0.30, 0.30, 0.128e3,
    0.01e9/100, 0.01e9/100, 0.01e9/100, 0.01e9/100, 0.01e9/100, 0.01e9/100, 2e5/100, 2e5/100, 2e5/100)
mat = [uni, double, gelcoat, nexus, balsa]

idx = [3, 4, 2]  # material index
t = [0.000381, 0.00051, 18*0.00053]
theta = [0, 0, 20]*pi/180
layup1 = Layer.(mat[idx], t, theta)
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

segments = [layup1, layup2, layup3, layup4]
xbreak = [0.0, 0.0041, 0.1147, 0.5366, 1.0]

idx = [1, 5, 1]
t = [38*0.00053, 0.003125, 38*0.00053]
theta = [0, 0, 0]*pi/180
web = Layer.(mat[idx], t, theta)

webs = [web, web]
webloc = [0.15, 0.5]
nodes, elements = afmesh(xaf, yaf, chord, twist, paxis, xbreak, webloc, segments, webs)

cache = initialize_cache(nodes, elements)
S, sc, tc = compliance_matrix(nodes, elements; cache)

F = [0.0; 0; 0]
M = [2e4; 1e4; 0]
strain_b, stress_b, strain_p, stress_p = strain_recovery(F, M, nodes, elements, cache)

failure = tsai_wu(stress_p, elements)

include("./pipe.jl")
nodes, elements = composite_pipe()

cache = initialize_cache(nodes, elements)
S, sc, tc = compliance_matrix(nodes, elements; cache)

F = [0.0; 0; 0]
M = [-1000.0; 0; 0]
strain_b, stress_b, strain_p, stress_p = strain_recovery(F, M, nodes, elements, cache)

idx = 481:500  # elements at x = 0 from y = 0.3 -> 0.5
n = length(idx)
yvec = zeros(n)
s11 = zeros(n)
s22 = zeros(n)
for i = 1:n
    _, _, yvec[i] = GXBeam.area_and_centroid_of_element(nodes[elements[idx[i]].nodenum])
    s11[i] = stress_b[3, idx[i]]
    s22[i] = stress_b[1, idx[i]]
end

data1 = [
-1.0408340855860843e-17  -0.2233363719234286
0.02504493708807669  -0.24156791248860665
0.049970041941282226  -0.25979945305378427
0.07501497902935894  -0.2734731084776677
0.10005991611743567  -0.29170464904284577
0.10005991611743564  -4.416590701914313
0.12510485320551235  -4.690063810391979
0.15002995805871785  -4.96809480401094
0.17507489514679453  -5.241567912488606
0.2001198322348712  -5.519598906107569
]

data2 = [0.00011487650775416497  0.004832104832104944
0.02504307869040781  0.06871416871416883
0.04997128087306147  0.13226863226863234
0.07501435956346927  0.19582309582309587
0.09994256174612295  0.2597051597051597
0.10005743825387711  -0.10196560196560189
0.12498564043653071  -0.10491400491400493
0.1500287191269386  -0.10786240786240764
0.17484204480183801  -0.11146601146601148
0.19999999999999993  -0.11408681408681415
]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

