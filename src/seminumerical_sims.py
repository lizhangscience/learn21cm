import numpy as np
import c2raytools as c2t
import usefuls

def sem_num(dens_cube, sourcelist, Nion=30, Nrec=0, Rmfp=10, boxsize=None):
	"""
	@Majumdar et al. (2014)
	"""
	if boxsize is None: boxsize = c2t.conv.LB
	Mhalos   = Mgrid_2_Msolar(sourcelist[:,3])
	xx,yy,zz = (sourcelist[:,:3].T-1).astype(dtype=np.int)
	N_h = np.zeros(dens_cube.shape)
	N_h[xx,yy,zz] = Nion*Mhalos*c2t.const.OmegaB/c2t.const.Omega0/c2t.const.m_p/c2t.const.solar_masses_per_gram
	nn  = dens_cube.shape[0]
	n_h = N_h*(nn/boxsize)**3
	n_H = c2t.const.Mpc**3*dens_cube/c2t.const.m_p
	G_mfp = Rmfp*nn/boxsize
	Rs = np.arange(G_mfp)+1.
	xf = np.zeros(dens_cube.shape)
	for i in xrange(Rs.size):
		ra     = Rs[i]
		kernel = spherical_tophat_kernel(nn, ra)
		nh_    = c2t.smooth_with_kernel(n_h, kernel)
		nH_    = c2t.smooth_with_kernel(n_H, kernel)
		xf[nh_>=nH_*(1+Nrec)] = 1
		print (i+1)*100/Rs.size, "% completed"
	xf[nh_<nH_*(1+Nrec)] = (nh_/nH_)[nh_<nH_*(1+Nrec)]
	return xf

def spherical_tophat_kernel(size, rad):
	kernel = np.zeros((size,size,size))
	center = np.array(kernel.shape)/2
	kernel = usefuls.put_sphere(kernel, center, rad, label=1)
	kernel /= np.sum(kernel)
	return kernel


	

