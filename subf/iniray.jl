function iniray(nr::Array{Int64,1},ns::Int64,nx::Int64,nz::Int64,
  dx::Float64,x0::Float64,rx::Array{Float64,2},sx::Array{Float64,1},
  slw::Array{Float64,1})

  ilflag = 2; # left path source
  irflag = 2; # left path source

  # find the global source/receiver boundary, the very left receiver
  # and its very left source; the very right receiver and its very
  # right source

  #  before input into this program, all sources and receivers are
  #  already sorted from left to right as a sequence

  if sx[ns] > sx[1]

    xmin = sx[1];
    ibs1 = 1;
    if rx[1,1] < rx[nr[1],1]
      ibr1 = nr[1];
    end
    if rx[1,1] > rx[nr[1],1]
      ibr1 = 1 ;
    end
    for is in 1:ns
      for ir in 1:nr[is]
        if rx[ir,is] <= xmin
          xmin = rx[ir,is];
          ibr1 = ir;
          ibs1 = is;
          ilflag = 1;    # left path receiver
        end
      end
    end

    xmax = sx[ns];
    ibs2 = ns;
    if rx[1,ns] < rx[nr[ns],ns]
      ibr2 = 1;
    end
    if rx[1,ns] > rx[nr[ns],ns]
      ibr2 = nr[ns];
    end
    for is in ns:-1:1
      for ir in 1:nr[is]
        if rx[ir,is] >= xmax
          xmax = rx[ir,is];
          ibr2 = ir;
          ibr2 = is;
          irflag = 1;     # right path receiver
        end
      end
    end

  else

    xmin = sx[ns];
    ibs1 = ns;
    if rx[1,ns] < rx[nr[ns],ns]
      ibr1 = nr[ns];
    end
    if rx[1,ns] > rx[nr[ns],ns]
      ibr1 = 1;
    end
    for is in ns:-1:1
      for ir in 1:nr[is]
        if rx[ir,is] <= xmin
          xmin = rx[ir,is]
          ibr1 = ir;
          ibs1 = is;
          ilflag = 1;
        end
      end
    end

    xmax = sx[1];
    ibs2 = 1;
    if rx[1,1] < rx[nr[1],1]
      ibr2 = 1;
    end
    if rx[1,1] > rx[nr[1],1]
      ibr2 = nr[1];
    end
    for is in 1:ns
      for ir in 1: nr[is]
        if rx[ir,is] >= xmax
          xmax = rx[ir,is];
          ibr2 = ir;
          ibs2 = is;
          irflag = 1;
        end
      end
    end

  end

  xx1 = 1.0e18;
  xx2 = -1.0e18;
  ixr1 = zeros(Float64,ns);
  ixr2 =zeros(Float64,ns);

  for is in 1:ns
#    xx1 = 1.0e18;            # !! notice here, we'll see the effects in the future
#    xx2 = -1.0e18;
    xx1 = min(xx1,sx[is]);
    xx2 = max(xx2,sx[is]);
    for ir in 1:nr[is]
      xx1 = min(xx1,rx[ir,is]);
      xx2 = max(xx2,rx[ir,is]);
    end
    ixr1[is] = Int64(floor((xx1 - x0) / dx + 1));  # source/receiver bound, grid, min
    ixr2[is] = Int64(floor((xx2 - x0) / dx + 1));  # source/receiver bound, grid, max
  end


  # remove air velocity if any for refraction case

  air = 340.0;
  for i in 1:nx
    for j = 1:nz
      ii = i + (j - 1) * nx;
      if slw[ii] < air
        for k in 1: j-1
          kk = i + (k - 1) * nx;
          slw[kk] = slw[ii];
        end
        break;
      end
    end
  end

  return ilflag,irflag,ibs1,ibs2,ibr1,ibr2,ixr1,ixr2,slw
end

# ibs1,ibs2,ibr1,ibr2, the source/receiver which is corresponding to
# the boundary, sequence is/ir

# ixr1,ixr2, the grid which is corresponding to the boundary, for
# each source
