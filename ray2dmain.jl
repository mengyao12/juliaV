using ray2dmod,PyPlot

parmpath = "/Users/admin/tomo2Dcode/juliaV/data/inparm.txt";
ttpth,mdlpth,outpth,maxitr,maxV,minV,tau,hvs = readparm(parmpath);

ns,nrall,mr,sx,sz,nr,rx,rz,dat = readgeom(ttpth);
nx,nz,dx,x0,z0,nm,slw = readmdl(mdlpth);

msg,ixff,ixrr,izup,izdn = checkin(nx,nz,x0,z0,dx,nr,ns,rx,rz,sx,sz);
msg,itop,ztop = mktop(nx,nz,x0,z0,dx,nr,ns,sx,sz,rx,rz,mr,ixff,ixrr,izup,izdn)
