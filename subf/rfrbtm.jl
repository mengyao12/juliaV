function rfrbtm(nx::Int64,nz::Int64,nm::Int64,slw::Array{Float64,1},
  izdn::Int64)

  # find the highest velocity vertically

  ibt = slwbtm(nx,nz,slw,izdn,nm);
  ix1 = 1;
  ix2 = 2;
  while ix1 < ix2
    apri = ibt[ix1+1] - ibt[ix1];
    for i in ix1+1:nx
      grd1 = (ibt[i] - ibt[ix1])/(i - ix1);
      if grd1 >= apri
        apri = grd1;
        ix2 = i;
      end
    end

    grd = (ibt[ix2] - ibt[ix1])/(ix2 - ix1);
    ibr = ibt[ix1];
    for k in ix1:ix2-1
      ibt[k] = ibr + 2 + Int64(floor(grd * (k - ix1)));
    end

    ix1 = ix2;
  end

  # ensure bottom safe

  for i in 1:nx
    if ibt[i] >= nz
      ibt[i] = nz - 1;
    end
  end

  return ibt
end
