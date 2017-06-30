using ray2dmod,PyPlot

parmpath = "/Users/admin/tomo2Dcode/juliaV/data/inparm.txt";

# set parameters
ms = 780;
mr = 1200;
mx = 3500;
mz = 400;
msmax = 80000;
mm = mx * mz;
imax = 4 * (mx + mz);
mxz = 1.2 * (mx + 2 * mz);

# read the input parameters
ttpth,mdlpth,outpth,maxitr,maxV,minV,tau,hvs = readparm(parmpath);

# read the geometry file
ns,nrall,mr,sx,sz,nr,rx,rz,dat = readgeom(ttpth);

# read the velocity model
nx,nz,dx,x0,z0,nm,slw = readmdl(mdlpth);

nxz = Int64(floor(1.2 * (nx + 2 * nz)));

# check the geometry boundary
msg,ixff,ixrr,izup,izdn = checkin(nx,nz,x0,z0,dx,nr,ns,rx,rz,sx,sz);

# build the topography according to the geometry
msg,itop,ztop = mktop(nx,nz,x0,z0,dx,nr,ns,sx,sz,rx,rz,mr,ixff,ixrr,izup,izdn);

# set the constant number for raytracing
sqr2,sqr5,sqr10,sqr13,sqr17,a4,c = inicon(dx);

# find the source/receive boundary, remove the air velocity
ilflag,irflag,ibs1,ibs2,ibr1,ibr2,ixr1,ixr2,slw = iniray(nr,ns,nx,nz,dx,x0,rx,
sx,slw);

# 2-D raytracing based on graph method
tt,ttime = graphray(nx,nz,nm,slw,izdn,itop,ixr1,ixr2,ns,dx,imax,x0,z0,sx,sz,
nxz,sqr2,sqr5,sqr10,sqr13,sqr17,a4,c,msmax,nr,rx,rz,mr);

#=
for ii in 1:nx*nz
  if tt[ii] == 1.0e18
    tt[ii] = 0.0;
  end
end
ttm = reshape(tt,nx,nz);
fig = imshow(ttm');
colorbar(fig)
=#
