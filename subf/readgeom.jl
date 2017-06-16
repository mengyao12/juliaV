function readgeom(ttpth::AbstractString)
  f = open(ttpth);
  ttm = readdlm(ttpth);
  ns = ttm[1,1];
  iunit=ttm[1,2];
  nr = zeros(Int64,ns);
  nrall = 0;
  mr = 0;
  cnt = 0;
  sx = zeros(Float64,ns);
  sz = zeros(Float64,ns);

  for i in 1:ns
    cnt = nrall + 2*i;
    nr[i] = ttm[cnt,2];
    nrall = nrall + nr[i];
    if nr[i] > mr
       mr = nr[i];
    end
    sx[i] = ttm[cnt+1,1];
    sz[i] = ttm[cnt+1,2];
  end

  rx = zeros(Float64,mr,ns);
  rz = zeros(Float64,mr,ns);
  dat = zeros(Float64,mr,ns);
  cnt = 0;
  nrall = 0;

  for i in 1:ns
    cnt = nrall + 2*i;
    nr[i] = ttm[cnt,2];
    nrall = nrall + nr[i];
    rx[1:nr[i],i] = ttm[cnt+2:cnt+2+nr[i]-1,2];
    rz[1:nr[i],i] = ttm[cnt+2:cnt+2+nr[i]-1,3];
    dat[1:nr[i],i] = ttm[cnt+2:cnt+2+nr[i]-1,4];
  end
  close(f)

  return ns,nrall,mr,sx,sz,nr,rx,rz,dat
end
