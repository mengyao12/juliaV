function timsor(nx::Int64,x0::Float64,z0::Float64,dx::Float64,is::Int64,
  sx::Array{Float64,1},sz::Array{Float64,1},nmele::Array{Int64,1},
  idele::Array{Int64,2},nm::Int64,ipath::Array{Int64,1},
  slw::Array{Float64,1},idone::Array{Int64,1},tt::Array{Float64,1})

  isx = Int64(floor((sx[is] - x0) / dx + 1));
  isz = Int64(floor((sz[is] - z0) / dx + 1));

  iss0 = isx + (isz - 1) * nx;
  xxs = (isx - 1) * dx + x0;
  zzs = (isz - 1) * dx + z0;
  ndid = 0;
  i1st = 1;
  nmele[1] = 0;
  i0 = Int64;
  j0 = Int64;


  if xxs == sx[is] && zzs == sz[is]   # right on a grid
    i = isx;
    j = isz;
    iss = i + (j - 1) * nx;
    tt[iss]=0.0;
    ipath[iss] = 0;
  else                               # not on a grid
    for i in isx+1:-1:isx
      for j in isz+1:-1:isz
        iss = i + (j - 1) * nx;
        if idone[iss] != 1
          tt[iss] = sqrt(((i - 1) * dx + x0 - sx[is])^2 +
          ((j - 1) * dx + z0 - sz[is])^2) * slw[iss0]
          ipath[iss] = 0;
          idone[iss] = 1;
          ndid = ndid + 1;
          nmele[1] = nmele[1] + 1;
          idele[nmele[1],1] = iss;
          i0 = i;
          j0 = j;
        end
      end
    end
    ndid = ndid - 1;
    i = i0;
    j = j0;
    iss = i + (j - 1) * nx;
  end

  return nmele,idele,ipath,idone,tt,i1st,ndid,iss,i,j
end
