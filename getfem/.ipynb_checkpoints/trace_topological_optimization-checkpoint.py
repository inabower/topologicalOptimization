import getfem as gf
import numpy as np

# parameters
NX=30;
ls_degree = 1;
alpha = 1;
beta = 1;
rayon_trous = 0.2;

m = gf.Mesh('cartesian', np.arange(-.5,.5+1./NX,1./NX), np.arange(-.5,.5+1./NX,1./NX))
mf_basic = gf.MeshFem(m, 1);
mf_basic.set_fem(gf.Fem('FEM_QK(2,2)'))
mls = gf.MeshLevelSet(m)
ls  = gf.LevelSet(m, ls_degree)
mls.add(ls)
mf_ls = ls.mf()

P = mf_ls.basic_dof_nodes()
x = P[0]
y = P[1]
ULS = 1000*np.ones(x.shape)

#% Loop on the topological optimization
#while True:
for col in range(1000):
    ls.set_values(ULS)
    mls.adapt()
    mim_bound = gf.MeshIm('levelset',mls,'boundary', gf.Integ('IM_TRIANGLE(6)'))
    mim = gf.MeshIm('levelset',mls,'outside', gf.Integ('IM_TRIANGLE(6)'))
    mim.set_integ(4)

    mf_mult = gf.MeshFem(m)
    mf_mult.set_fem(gf.Fem('FEM_QK(2,1)'))
    
    M = gf.asm_mass_matrix(mim, mf_basic)
    D = np.abs(M.diag().T[0])
    ind = np.argwhere(D>1e-8)
    mf = gf.MeshFem('partial', mf_basic, ind)

    S = gf.asm('volumic','V()+=comp()',mim)
    print('remaining surface :',S)

    # % Problem definition (Laplace(u) + u = f)
    md = gf.Model('real')
    md.add_fem_variable('u', mf)
    md.add_Laplacian_brick(mim, 'u')
    md.add_fem_data('VolumicData', mf_basic)
    md.add_source_term_brick(mim, 'u', 'VolumicData')
    md.add_initialized_data('rho', [1.0])
    md.add_mass_brick(mim, 'u', 'rho')
    md.add_multiplier('mult_dir', mf_mult, 'u')
    # % To be completely robust, a stabilization should be used on the Dirichlet
    # % boundary to ensure the inf-sup condition (Nitsche or Barbosa-Hughes)
    md.add_Dirichlet_condition_with_multipliers(mim_bound, 'u', 'mult_dir', -1)

    print('x shape :',x.shape)
    # % Solving the direct problem.
    U0 = mf_basic.eval('0.4*(3.*np.sin(np.pi*(x+y)) + ((x-0.5)**10 + (y-0.5)**10 + (x+0.5)**10 + (y+0.5)**10))', globals(), locals())
    U0 = np.array([U0])
    md.set_variable('VolumicData', U0)

    # gf_model_get(md, 'solve');
    md.solve()
    U = md.variable("u")

    """
    subplot(2,1,1);
    gf_plot(mf, U);
    hold on;
    [h1,h2]=gf_plot(mf_ls, get(ls,'values'), 'contour', 0,'pcolor','off');
    set(h2{1},'LineWidth',2);
    set(h2{1},'Color','green');
    colorbar;
    title('u');
    hold off;
    """
    
    # % Solving the adjoint problem.
    UBASIC = gf.compute_interpolate_on(mf, U, mf_basic)
    F = 2*UBASIC
    md.set_variable('VolumicData', F)
    md.solve()
    W = md.variable("u")

    # % Computation of the topological gradient
    mf_g = gf.MeshFem(m, 1)
    mf_g.set_fem(gf.Fem('FEM_PRODUCT(FEM_PK_DISCONTINUOUS(1,2),FEM_PK_DISCONTINUOUS(1,2))'))
    DU = gf.compute(mf, U, 'gradient', mf_g);
    DW = gf.compute(mf, W, 'gradient', mf_g);
    nbdof = mf_g.nbdof()
    DU = np.reshape(DU, [2, nbdof])
    DW = np.reshape(DW, [2, nbdof])
    
    UU = gf.compute_interpolate_on(mf, U, mf_g)
    UU0 = gf.compute_interpolate_on(mf_basic, U0, mf_g)
    LS = gf.compute_interpolate_on(mf_ls, ULS, mf_g)

    G = (-4*np.pi*( alpha*(DU[0]**2 + DU[1]**2 + DU[0]*DW[0] + DU[1]*DW[1]) + beta*(UU-UU0)**2)) * (np.sign(LS)+1.)/2;

    val = G.min()
    i   = G.argmin()
    point = mf_g.basic_dof_nodes(i)
    print('val =',val)
    
    if val >= -12:
        break

    R = -(val+7) / 200;
    xc = point[0][0]
    yc = point[1][0]
    P = mf_ls.basic_dof_nodes()
    x = P[0]
    y = P[1]
    ULS = np.minimum(ULS, ((x - xc)**2 + (y - yc)**2) - R**2);
    mf.export_to_vtk(f'fig/U_{col:03d}.vtk', U, 'U')
    mf_g.export_to_vtk(f'fig/G_{col:03d}.vtk', G, 'G')
    mf_g.export_to_vtk(f'fig/LS_{col:03d}.vtk', LS, 'LS')
    mf_ls.export_to_vtk(f'fig/ULS_{col:03d}.vtk', ULS, 'ULS')