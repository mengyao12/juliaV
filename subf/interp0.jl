function interp0(iirc::Int64,iirp::Int64,nod::Array{Int64,1},dx::Float64,
  nx::Int64,nz::Int64)

  # tell nodes and ray segments between "iirc" and "iirp"
  # here, nod represent the corresponding "source" (iirp) of current
  # "receiver" (iirc)
  # inod represent the number of nodes that a ray pass through

  i0,j0 = denum(iirp,nx,nz);
  i1,j1 = denum(iirc,nx,nz);
  mm = (i1 - i0)^2 + (j1 - j0)^2;

  if mm == 1
    inod = 1;
    nod[1] = iirp;
    if i1 < i0 || j1 < j0
      nod[1] = iirc;     # plot the grid to see which ray segment is prefered
    end
      return inod,nod
  end

  if mm == 2
    inod = 1;
    if i0 < i1 && j0 < j1
      nod[1] = iirp;
    end
    if i0 > i1 && j0 < j1
      nod[1] = iirp - 1;
    end
    if i0 < i1 && j0 > j1
      nod[1] = iirc - 1;
    end
    if i0 > i1 && j0 > j1
      nod[1] = iirc;
    end
      return inod,nod
  end

    if mm == 5
      inod = 2;
    end
    if mm == 10
      inod = 3;
    end
    if mm == 13
      inod = 4;
    end
    if mm == 17
      inod = 4;
    end
    if mm == 25
      inod =6;
    end

    if iirc > iirp
      ii1 = i1;
      jj1 = j1;
      ii0 = i0;
      jj0 = j0;
    else
      ii1 = i0;
      jj1 = j0;
      ii0 = i1;
      jj0 = j1;
    end

    k = 1;
    if ii1 > ii0 && jj1 > jj0
      ip = 0;
    @label path1
      nod[k] = ii0 + (jj0 - 1) * nx;
      if ii0 == ii1-1 && jj0 == jj1-1
        @goto pathend
      end
      if ip == 0
        @goto path2
      end
      if ip == 1
        @goto path3
      end
    @label path2
      if jj0+1 < jj1
        jj0 = jj0 + 1;
        k = k + 1;
      end
      ip = 1;
      @goto path1
    @label path3
      if ii0+1 < ii1
        ii0 = ii0 + 1;
        k = k + 1;
      end
      ip = 0;
      @goto path1
    end

    if ii1 < ii0 && jj1 > jj0
      ip = 0;
      ii0 = ii0 - 1;
    @label path4
      nod[k] = ii0 + (jj0-1) * nx;
      if ii0 == ii1 && jj0+1 == jj1
        @goto pathend
      end
      if ip == 0
        @goto path5
      end
      if ip == 1
        @goto path6
      end
    @label path5
      if jj0 + 1 < jj1
        jj0 = jj0 + 1;
        k = k + 1;
      end
      ip = 1;
      @goto path4
    @label path6
      if ii0-1 >= ii1
        ii0 = ii0 - 1;
        k = k + 1;
      end
      ip = 0;
      @goto path4
    end

    @label pathend

    return inod,nod
  end
