function slwminmax(nx::Int64,nm::Int64,slw::Array{Float64,1},
  itop::Array{Int64,1},ibt::Array{Int64,1})

  icnt = 0;
  wrk1 = zeros(Float64,nm);
  for i in 1:nx
    for j in itop[i]:ibt[i]
      iss = i + (j - 1) * nx;
      icnt = icnt + 1;
      wrk1[icnt] = slw[iss];
    end
  end

  slwmin = minimum(wrk1[1:icnt]);
  slwmax = maximum(wrk1[1:icnt]);

  return slwmin,slwmax
end
