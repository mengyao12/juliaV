function line(x1::Float64,x2::Float64,z1::Float64,z2::Float64,xx::Float64)

  grd = (z2 - z1) / (x2 - x1);
  zz = z2 + grd * (xx - x2);

  return zz
end
