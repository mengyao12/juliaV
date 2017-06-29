function denum(iss::Int64,i::Int64,j::Int64,nx::Int64,nz::Int64)

  j = Int64(floor(iss / nx + 1));
  i = Int64(floor(iss - (j - 1) * nx));
  if i == 0
    i = nx;
    j = j - 1;
  end

  return i,j
end
