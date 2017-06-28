function washtomo(is::Int64,ns::Int64,nx::Int64,nm::Int64,
  ixr1::Array{Int64,1},ixr2::Array{Int64,1},itop::Array{Int64,1},
  ibt::Array{Int64,1},idone::Array{Int64,1},tt::Array{Float64,1})

  for i in ixr1[is]-1:ixr2[is]+1
    for j in itop[i]:ibt[i]
      ii = i + (j-1) * nx;
      idone[ii] = 0;
    end
  end

  for i in 1:nm
    tt[i] = 1.0e18;
  end

  return idone,tt
end
