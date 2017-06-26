function slwbtm(nx::Int64,nz::Int64,slw::Array{Float64,1},izdn::Int64,nm::Int64)

  # find the highest velocity vertically

  ibt = zeros(Int64,nx);
  for i in 1:nx
    s0 = slw[i+izdn*nx];
    ibt[i] = izdn + 1;
    for j in izdn+2:nz-1
      ii = i + (j - 1) * nx;
      if slw[ii] < s0
        ibt[i] = j;
        s0 = slw[ii];
      end
    end
  end
      
  return ibt
end
