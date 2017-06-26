function readmdl(mdlpth::AbstractString)
  f = open(mdlpth);
  mdlm = readdlm(mdlpth);
  nx = mdlm[1,1];
  nz = mdlm[1,2];
  dx = mdlm[1,3];
  x0 = mdlm[1,4];
  z0 = mdlm[1,5];
  iunit = mdlm[1,6];
  nm = nx * nz

  slw = zeros(Float64,nx*nz);
  mx,mz = size(mdlm);
  cnt=0;
  for i in 2:mx
    for j in 1:mz
    if typeof(mdlm[i,j]) != SubString{String}
      cnt = cnt+1;
      slw[cnt] = 1.0/mdlm[i,j];
    end
    end
  end
  close(f)

return nx,nz,dx,x0,z0,nm,slw

end
