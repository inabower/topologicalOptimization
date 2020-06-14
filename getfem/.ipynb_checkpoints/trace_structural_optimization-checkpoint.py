import getfem as gf
import numpy as np
from scipy.sparse.linalg import cg
import time

TEST_CASE = 1
if TEST_CASE == 0:
    N = 2;
    initial_holes = 1;
elif TEST_CASE == 1:
    N = 2;
    initial_holes = 0;
elif TEST_CASE == 2:
    N = 3;
    initial_holes = 1;
elif TEST_CASE == 3:
    N = 3;
    initial_holes = 0;
    
k = 1;           #% Degree of the finite element method for u
λ = 1;           #% Lame coefficient
μ = 1;           #% Lame coefficient

if N == 2:
    NY = 80             # % Number of elements in y direction
    level_set_rate = 0.4 / NY
    reinitialisation_time = 0.005;
    threshold_shape = 0.90;
    if TEST_CASE == 1:
        threshold_topo  = 1.3;
    else:
        threshold_topo  = 0;
    penalty_param = 1e-6;
    nbiter = 400;
    NBDRAW = 20;            #% Draw solution each NBDRAW iterations
else:
    NY = 30
    level_set_rate = 0.025 / NY;
    reinitialisation_time = 0.0035;
    threshold_shape = 28;
    if TEST_CASE == 3:
        threshold_topo = 30;
    else:
        threshold_topo = 0;
    penalty_param = 1e-6;
    nbiter = 600;
    NBDRAW = 5;            #% Draw solution each NBDRAW iterations
    
hole_radius = max(0.03, 2./NY);  #% Hole radius for topological optimization
ls_degree = 1;         #% Degree of the level-set. Should be one for the moment.
DEBUG = False
if DEBUG:
    NG = 3;
else:
    NG = 2;
    
#% Mesh definition
#% m=gfMesh('cartesian', -1:(1/NY):1, -.5:(1/NY):.5);
if (N == 2):
    m=gf.Mesh('regular simplices', np.linspace(-1.,1.,NY), np.linspace(-.5,.5,NY))
else:
    m=gf.Mesh('regular simplices', np.linspace(-1.,1.,NY), np.linspace(-.5,.5,NY), np.linspace(-.5,.5,NY));
pts =  m.pts()

#% Find the boundary GammaD and GammaN
pidleft = np.compress((abs(pts[0,:]+1.0)<1e-7), list(range(0, m.nbpts())))
fidleft = m.faces_from_pid(pidleft)
normals = m.normal_of_faces(fidleft)
fidleft = fidleft[:, abs(normals[0,:]+1.0)<1e-3]
GAMMAD = 2;
m.set_region(GAMMAD, fidleft)

pidright = np.compress((abs(pts[0,:]-1.0)<1e-7), list(range(0, m.nbpts())))
fidright = m.faces_from_pid(pidright)
normals = m.normal_of_faces(fidright)
fidright = fidright[:, abs(normals[0,:]-1.0)<1e-3]
GAMMAN = 3;
m.set_region(GAMMAN, fidright)

#% Definition of the finite element methods
ls  = gf.LevelSet(m, ls_degree)
mls = gf.MeshLevelSet(m)
mls.add(ls)
mf_ls = ls.mf()
if N==2:
    mimls = gf.MeshIm(m, gf.Integ('IM_TRIANGLE(4)'))
else:
    mimls = gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(5)'))

mf_basic = gf.MeshFem(m, N)
mf_basic.set_fem(gf.Fem(f'FEM_PK({N},{k})'))

mf_g = gf.MeshFem(m, 1)
mf_g.set_fem(gf.Fem(f'FEM_PK_DISCONTINUOUS({N},{k-1})'))
        
mf_cont = gf.MeshFem(m, 1)
mf_cont.set_fem(gf.Fem(f'FEM_PK({N},{ls_degree})'))

print(f'There is {mf_basic.nbdof()} elasticity dofs')

Mcont = gf.asm_mass_matrix(mimls, mf_cont)
RMcont = np.linalg.cholesky(Mcont.full())
Mcontls = gf.asm_mass_matrix(mimls, mf_ls)
RMcontls = np.linalg.cholesky(Mcontls.full())

#% Definition of the initial level-set
if (initial_holes):
    if (N == 2):
        ULS = mf_basic.eval('(-0.6-np.sin(np.pi*4*x)*np.cos(np.pi*4*y))/(4*np.pi)', globals(), locals())
    else:
        ULS = mf_basic.eval('(-0.6-np.sin(np.pi*4*x)*np.cos(np.pi*4*y)*np.cos(np.pi*4*z))/(4*np.pi)', globals(), locals())
else:
    ULS = mf_basic.eval('x - 2', globals(), locals())
    
P = mf_ls.basic_dof_nodes()
#% Force on the right part (Neumann condition)
if (N == 2):
    F = mf_basic.eval('[0, -1.0*(abs(y) < 0.05)]', globals(), locals())
else:
    F = mf_basic.eval('[0,0,-20*(np.abs(y) < 0.05)*(np.abs(z) < 0.05)]', globals(), locals())
    
ls.set_values(ULS)
print('Adapting the mesh')
mls.adapt()
if N == 2:
    mim = gf.MeshIm('levelset',mls,'inside', gf.Integ('IM_TRIANGLE(6)'))
else:
    mim = gf.MeshIm('levelset',mls,'inside', gf.Integ('IM_TETRAHEDRON(6)'))
print('Mesh adapted')
mim.set_integ(4)
print('Integration methods adapted')
mf = gf.MeshFem('partial', mf_basic, range(int(mf_basic.nbdof()/2)))

md = gf.Model('real')
md.add_fem_variable('u', mf)
md.add_initialized_data('mu', [μ])
md.add_initialized_data('lambda', [λ])
md.add_isotropic_linearized_elasticity_brick(mim, "u", "lambda", "mu")
md.add_Dirichlet_condition_with_multipliers(mim, 'u', 1, GAMMAD)
md.add_Dirichlet_condition_with_penalization(mim, 'u', 1e9, GAMMAD)
md.add_initialized_data('penalty_param', [penalty_param])
md.add_mass_brick(mim, 'u', 'penalty_param')
md.add_initialized_fem_data("Force", mf_basic, F)
md.add_source_term_brick(mim, 'u', 'Force', GAMMAN)

for niter in range(nbiter):
    start = time.time()
    if niter > 0:
        ls.set_values(ULS)
        print('Adapting the mesh')
        mls.adapt()
        print('Mesh adapted')
        mim.adapt()
        print('Integration methods adapted')
    Mcont = gf.asm_mass_matrix(mimls, mf_cont)
    #D = np.abs(np.diag(np.diag(Mcont.full())))
    #ind = np.argwhere(D > ((1/NY)**N)*1e-7)
    D = np.abs(Mcont.diag().T[0])
    ind = np.argwhere(D > ((1/NY)**N)*1e-7)
    mf = gf.MeshFem('partial', mf_basic, ind)
    
    #% Solving the direct problem
    print('solving the direct problem')
    
    md.solve('max_res', 1e-7)
    U = md.variable("u")
    nbd = mf_ls.nbdof()
    
    #% Computation of indicators (computation of K could be avoided)
    K = gf.asm('linear elasticity', mim, mf, mf_ls, λ*np.ones((1, nbd)), μ*np.ones((1, nbd)));
    print(f'Elastic energy at iteration {niter}: {np.dot(U,(K*U))}')
    S = gf.asm('volumic','V()+=comp()',mim);
    if (N == 2):
        print(f'Remaining surface of material: {S}')
    else:
        print(f'Remaining volume of material: {S}')
        
    DU = gf.compute(mf, U, 'gradient', mf_g);
    EPSU = DU + DU[[1,0],:]
    
    #% Computation of the shape derivative
    if (N == 2):
        GF1 = (DU[0][0] + DU[1][1])**2*λ + 2*μ*np.sum(np.sum(EPSU, axis=0), axis=0);
    else:
        GF1 = (DU[0][0] + DU[1][1] + DU[2][2])**2*λ + 2*μ*np.sum(np.sum(EPSU, axis=0), axis=0);
    GF = GF1.reshape((1, len(GF1))) - threshold_shape;
    
    #% computation of the topological gradient
    if (N == 2):
        GT = -np.pi*( (λ+2*μ) / (2*μ*(λ+μ)) * (4*μ*GF1 + 2*(λ-μ)*(λ+μ)*(DU[0][0] + DU[1][1])**2));
    else:
        GT = -np.pi*( (λ+2*μ) / (μ*(9*λ+14*μ)) * (20*μ*GF1 + 2*(3*λ-2*μ)*(λ+μ)*(DU[0][0] + DU[1][1] + DU[2][2])**2));
    GT = GT.reshape((1, len(GT))) + threshold_topo;
    
    M = gf.asm_mass_matrix(mim, mf_g)
    D = np.abs(M.diag().T[0])
    maxD = np.max(D);
    ind = np.argwhere(D < maxD/40)
    mf = gf.MeshFem('partial', mf_basic, ind)
    GF[ind] = 0.0
    GF = np.clip(GF, None, 2*threshold_shape);
    ind = np.argwhere(D < maxD/1.2)
    GT[ind] = -20.0
    
    #% Drawing the gradients
    
    i   = GT.argmin()
    val = GT.min()
    print(f'Max value of the topological gradient: {val}')
    if (val > 0):
        point = mf_g.basic_dof_nodes(i)
        if (N == 2) :
            print(f'Making a new hole whose center is ({point[0]}, {point[1]})')
            ULS = np.max(ULS, (hole_radius**2 - (P[0,:] - point[0])**2 - (P[1,:] - point[1])**2)/(2*hole_radius))
        else:
            print(f'Making a new hole whose center is ({point[0]}, {point[1]}, {point[2]})')
            ULS = np.max(ULS, (hole_radius**2 - (P[0,:] - point[0])**2 - (P[1,:] - point[1])**2 - P[2,:] - point[2]**2)/(2*hole_radius))
    Fdisc = gf.asm('volumic source', mimls, mf_ls, mf_g, GF);
    #cgs(A,b,tol,maxit,M1,M2,x0)
    #scipy.sparse.linalg.cg(A, b, x0=None, tol=1e-05, maxiter=None, M=None, callback=None, atol=None)
    Vcont = cg(Mcontls.full(), Fdisc, tol=1e-8, maxiter=1000, M=RMcontls);
    ULS = ULS + Vcont[0] * level_set_rate
    
    #% Re-initialization of the level-set
    dt = reinitialisation_time
    NT = 10
    ddt = dt / NT
    ULS0 = ULS
    
    for t in np.arange(ddt,ddt,dt):
        DLS = gf.compute(mf_ls, ULS, 'gradient', mf_g);
        Fdisc = gf.asm('volumic source', mimls, mf_cont, mf_g, DLS);
        DLScont = cg(Mcont.full(), Fdisc, tol=1e-8, maxiter=1000, M=RMcont)
        NORMDLS = np.sqrt(np.sum(np.reshape(DLScont, [N, len(DLScont[0])/N])**2,1))+1e-12;
        SULS = np.sign(ULS0) / NORMDLS;
        
        if (N == 2):
            W = DLScont*np.reshape([SULS, SULS], [N*len(SULS[1]), 1]);
        else:
            W = DLScont*np.reshape([SULS, SULS, SULS], [N*len(SULS[1]), 1]);
            
        gf.compute(mf_ls, ULS, 'convect', mf_cont, W, ddt, 1, 'unchanged');
        ULS = ULS + ddt * np.sign(ULS0);
    t = time.time() - start
    print(f'this iteration took {int(t/60)} minutes');