function checkin(nx::Int64,nz::Int64,x0::Float64,z0::Float64,dx::Float64,
  nr::Array{Int64,1},ns::Int64,rx::Array{Float64,2},rz::Array{Float64,2},
  sx::Array{Float64,1},sz::Array{Float64,1})

  if ns <= 0
    msg = "Error: Number of sources of zero!";
    return msg
  end

  iis = 0;
  for is in 1:ns-1
    ddx = abs(sx[is+1] - sx[is]);
    ddz = abs(sz[is+1] - sz[is]);
    if ddz > ddx     # shots in hole
      iis = iis + 1;
    end
  end

  iir = 0;
  for ir in 1:nr[1]-1
    ddx = abs(rx[ir+1,1] - rx[ir,1]);
    ddz = abs(rz[ir+1,1] - rz[ir,1]);
    if ddz > ddx      # receivers in hole
      iir = iir + 1;
    end
  end

  if iis == ns-1 || iir == nr[1]-1
    msg = "This survey is not in refraction geometry!";
    return msg
  end

  ainf = 1.0e18;
  xff = ainf;      # left
  xrr = - ainf;    # right
  zup = ainf;      # up
  zdn = - ainf;    # down
  for is in 1:ns
    xff = min(xff,sx[is]);
    xrr = max(xrr,sx[is]);
    zup = min(zup,sz[is]);
    zdn = max(zdn,sz[is]);
    for ir in 1:nr[is]
      xff = min(xff,rx[ir,is]);
      xrr = max(xrr,rx[ir,is]);
      zup = min(zup,rz[ir,is]);
      zdn = max(zdn,rz[ir,is]);
    end
  end

  ixff = Int64(floor((xff - x0) / dx + 1));
  ixrr = Int64(floor((xrr - x0) /dx + 1));
  izup = Int64(floor((zup - z0) / dx + 1));
  izdn = Int64(floor((zdn - z0) /dx + 1));

  if ixff <= 1
    msg = "Geometry beyond model left edge!";
    return msg
  end

  if ixrr >= nx
    msg = "Geometry beyong model right edge!";
    return msg
  end

  if izup <= 1
    msg = "Geometry is above model top boundary!";
    return msg
  end

  if izdn >= nz
    msg = "Geometry is below model bottom boundary!";
    return msg
  end

  msg = "This is a correct refraction geometry survey.";
  return msg,ixff,ixrr,izup,izdn
end
