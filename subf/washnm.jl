function washnm(nmele::Array{Int64,1},imax::Int64)

  for i in 1:imax
    nmele[i] = 0;
  end

  return nmele
end
