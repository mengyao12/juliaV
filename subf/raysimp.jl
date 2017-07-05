function raysimp(ij0::Int64,i0::Int64,j0::Int64,ipath::Array{Int64,1},nx::Int64,
  nz::Int64)

  ij1 = ipath[ij0];
  if ij1 == 0
    @goto raysimpend
  end

@label raylp1
  i1,j1 = denum(ij1,nx,nz);
@label raylp2
  ij2 = ipath[ij1];
  if ij2 == 0
    @goto raysimpend
  end

  i2,j2 = denum(ij2,nx,nz);

  kk1 = (i1 - i0) * (j2 - j1);
  kk2 = (i2 - i1) * (j1 - j0);

  if kk1 == kk2
    ipath(ij0) = ij2;
    ij1 = ij2;
    i1 = i2;
    j1 = j2;
    @goto raylp2
  else
    ij0 = ij1;
    i0 = i1;
    j0 = j1;
    @goto raylp1
  end

  @label raysimpend

  return ipath,i0,j0
