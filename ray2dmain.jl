using ray2dmod,PyPlot

parmpath = "/Users/admin/tomo2Dcode/juliaV/data/inparm.txt";
ttpth,mdlpth,outpth,maxitr,maxV,minV,tau,hvs = readparm(parmpath);

ns,nrall,mr,sx,sz,nr,rx,rz,dat = readgeom(ttpth);
nx,nz,dx,x0,z0,nm,slw = readmdl(mdlpth);

# check the geometry boundary
msg,ixff,ixrr,izup,izdn = checkin(nx,nz,x0,z0,dx,nr,ns,rx,rz,sx,sz);

# build the topography according to the geometry
msg,itop,ztop = mktop(nx,nz,x0,z0,dx,nr,ns,sx,sz,rx,rz,mr,ixff,ixrr,izup,izdn);

# set the constant number for raytracing
sqr2,sqr5,sqr10,sqr13,sqr17,a4,c = inicon(dx);

# find the source/receive boundary, remove the air velocity
ilflag,irflag,ibs1,ibs2,ibr1,ibr2,ixr1,ixr2,slw = iniray(nr,ns,nx,nz,dx,x0,rx,sx,slw);

ibt = slwbtm(nx,nz,slw,izdn,nm);
