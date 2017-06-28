function update(tt0::Float64,nd0::Int64,tt::Array{Float64,1},
  ipath::Array{Int64,1},ndid::Int64,iss::Int64,tmin::Float64,
  nmele::Array{Int64,1},idele::Array{Int64,2})

  if tt0 < tt[nd0]
    tt[nd0] = tt0;
    ipath[nd0] = iss;
    ndid = ndid + 1;
    ii = Int64(floor((tt0 / tmin)+1));
    nmele[ii] = nmele[ii] + 1;
    idele[nmele[ii],ii] = nd0;
  end

  return tt,ipath,ndid,nmele,idele
end
