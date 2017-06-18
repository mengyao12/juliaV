function mktop(nx::Int64,nz::Int64,x0::Float64,z0::Float64,dx::Float64,
  nr::Array{Int64,1},ns::Int64,sx::Array{Float64,1},sz::Array{Float64,1},
  rx::Array{Float64,2},rz::Array{Float64,2},mr::Int64,ixff::Int64,
  ixrr::Int64,izup::Int64,izdn::Int64)

  ztop = zeros(Float64,nx);
  itop = zeros(Int64,nx);
  ainf = 1.0e18;
  for i in 1:nx
    ztop[i] = ainf;
    itop[i] = ainf;
  end

  # assign topography value at source and receiver locations, perfer
  # receiver depth if it overlaps with a source

  for is in 1:ns
    ii = Int64(floor((sx[is] - x0) / dx + 1));
    ztop[ii] = sz[is];
  end

  for is in 1:ns
    for ir in 1:nr[is]
      ii = Int64(floor((rx[ir,is] - x0) / dx +1));
      if rz[ir,is] < ztop[ii]
        ztop[ii] = rz[ir,is];
      end
    end
  end

  # interpolate from left to right

  ixx1 = ixff;
  for ii in ixx1+1:ixrr
    if ztop[ii] < ainf
      ixx2 = ii;
      x1 = Float64(ixx1);
      x2 = Float64(ixx2);
      z1 = ztop[ixx1];
      z2 = ztop[ixx2];
      for jj in ixx1+1:ixx2-1
        xx = Float64(jj);
        zz = line(x1,x2,z1,z2,xx);
        ztop[jj] = zz;
      end
    else
      continue
    end
      if ixx1 != ixrr
        ixx1 = ixx2;
      end
  end

  # left edge extrapolation

  for ii in 1:ixff-1
    ztop[ii] = ztop[ixff];
  end

  # right edge extrapolation

  for ii in ixrr+1: nx
    ztop[ii] = ztop[ixrr];
  end

  # check source/receiver depth, if any under topography

  for is in 1:ns
    isx = Int64(floor((sx[is] - x0) / dx +1));
    isz = Int64(floor((sz[is] - z0) / dx +1));
    if isz < itop[isx]
      itop[isx] = isz;
    end
      for ir in 1:nr[is]
        irx = Int64(floor((rx[ir,is] - x0) / dx + 1));
        irz = Int64(floor((rz[ir,is] - z0) / dx + 1));
        if irz < itop[irx]
          itop[irx] = irz;
        end
        if irz < itop[irx+1]
          itop[irx+1] = irz;
        end
      end
  end

  # check topography depth
  for i in 1:nx
    itop[i] =Int64(floor((ztop[i] - z0) / dx +1));
    if itop[i] <= 1
      msg = "Topography is shallower than model top!";
      return msg
    end
    if itop[i] >= nz
      msg = "Topography is deeper than model bottom!";
      return msg
    end
  end

  msg = "Topography is correctly built."

return msg,itop,ztop
end
