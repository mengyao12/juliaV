function timrec(is::Int64,nr::Array{Int64,1},nx::Int64,x0::Float64,
  z0::Float64,dx::Float64,rx::Array{Float64,2},rz::Array{Float64,2},
  tt::Array{Float64,1},ttime::Array{Float64,2})

  for ir in 1:nr[is]
    irx = Int64(floor((rx[ir,is] - x0) / dx + 1));
    irz = Int64(floor((rz[ir,is] - z0) / dx + 1));
    xxr = (irx - 1) * dx + x0;
    zzr = (irz - 1) * dx + z0;

    t1 = tt[irx + (irz - 1) * nx];
    t2 = tt[irx + 1 + (irz - 1) * nx];
    t3 = tt[irx + irz * nx];
    t4 = tt[irx + 1 + irz * nx];
    xx = rx[ir,is] - xxr;
    zz = rz[ir,is] - zzr;
    t = xx / dx;
    u = zz / dx;

    ttime[ir,is] = (1.0 - t)*(1.0 - u) * t1 + t * (1.0 - u) *t2 +
    t * u * t3 + (1.0 - t) * u * t4;
  end

  return ttime
end
