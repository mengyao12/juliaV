function graphray(nx::Int64,nz::Int64,nm::Int64,slw::Array{Float64,1},
  izdn::Int64,itop::Array{Int64,1},ixr1::Array{Int64,1},
  ixr2::Array{Int64,1},ns::Int64,dx::Float64,imax::Int64,
  x0::Float64,z0::Float64,sx::Array{Float64,1},sz::Array{Float64,1},
  nxz::Int64,sqr2::Float64,sqr5::Float64,sqr10::Float64,sqr13::Float64,
  sqr17::Float64,a4::Float64,c::Array{Float64,1},msmax::Int64,
  nr::Array{Int64,1},rx::Array{Float64,2},rz::Array{Float64,2},mr::Int64,
  iorder::Int64)

  # ensure the raytracing bottom area
  ibt = rfrbtm(nx,nz,nm,slw,izdn);

  # get the min and max slw of the velocity model
  slwmin,slwmax = slwminmax(nx,nm,slw,itop,ibt);

  tmin = dx * slwmin * 0.995;
  nmele = zeros(Int64,imax);
  idele = zeros(Int64,nxz,imax);
  ipath = zeros(Int64,nm);
  indcn = zeros(Int64,nm);
  idone = zeros(Int64,nm);
  tt = zeros(Float64,nm);
  nod = zeros(Int64,200);
  i1st = Int64;
  ndid = Int64;
  i = Int64;
  j = Int64;
  iss = Int64;
  tt0 = Float64;
  iirp = Int64;
  iirc = Int64;
  inod = Int64;
  nd = zeros(Int64,6);
  icnt = Int64;
  ip = Int64;
  ip = 1;
  iasen = zeros(Int64,msmax*ns);
  ttime = zeros(Float64,mr,ns)

  for ii in 1:nm
    idone[ii] = 1;
  end

  for is in 1:2   # !! source loop start
    println(is)

    # initialize the grid for each source
    idone,tt = washtomo(is,ns,nx,nm,ixr1,ixr2,itop,ibt,idone,tt);
    nmele = washnm(nmele,imax);
    indcn = washnm(indcn,nm);

    # calculate the tt of source
    nmele,idele,ipath,idone,tt,i1st,ndid,iss,i,j = timsor(nx,x0,z0,dx,is,sx,sz,
    nmele,idele,nm,ipath,slw,idone,tt);

    # order 2, 8 nodes
@label order
      indcn[iss] = 1;

      if j >= 2

        if i >= 2
          nod[1] = i - 1 + (j - 2) * nx;   # upper row 2 1st node
          if idone[nod[1]] != 1
            tt0 = sqr2 * slw[nod[1]] + tt[iss];
            tt,ipath,ndid,nmele,idele  = update(tt0,nod[1],tt,ipath,ndid,iss,tmin,
            nmele,idele);
          end
        end
        if i >= 2 && i <= nx-1
          nod[2] = i + 0 + (j - 2) * nx;  # upper row 2 2nd node
          if idone[nod[2]] != 1
            tt0 = dx * (min(slw[nod[1]],slw[nod[2]])) + tt[iss];
            tt,ipath,ndid,nmele,idele  = update(tt0,nod[2],tt,ipath,ndid,iss,tmin,
            nmele,idele);
          end
        end
        if i <= nx-1
          nod[3] = i + 1 + (j - 2) * nx;  # upper row 2 3rd node
          if idone[nod[3]] != 1
            tt0 = sqr2 * slw[nod[2]] + tt[iss];
            tt,ipath,ndid,nmele,idele  = update(tt0,nod[3],tt,ipath,ndid,iss,tmin,
            nmele,idele);
          end
        end

      end

      if i >= 2
        nod[4] = i - 1 + (j - 1) * nx;     # middle 1st node
        if idone[nod[4]] != 1
          tt0 = dx * min(slw[nod[1]],slw[nod[4]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[4],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end
      if i <= nx-1
        nod[5] = i + 1 + (j - 1) * nx;    # middle 2nd node
        if idone[nod[5]] != 1
          tt0 = dx * (min(slw[nod[2]],slw[iss])) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[5],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if j <= nz-1

        if i >= 2
          nod[6] = i - 1 + (j + 0) * nx;   # lower row 2 1st node
          if idone[nod[6]] != 1
            tt0 = sqr2 * slw[nod[4]] + tt[iss];
            tt,ipath,ndid,nmele,idele  = update(tt0,nod[6],tt,ipath,ndid,iss,tmin,
            nmele,idele);
          end
        end
        nod[7] = i + 0 + (j + 0) * nx;    # lower row 2 2nd node
        if idone[nod[7]] != 1
          tt0 = dx * (min(slw[nod[4]],slw[iss])) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[7],tt,ipath,ndid,iss,tmin,
        nmele,idele);
        end
        if i <= nx-1
          nod[8] = i + 1 + (j + 0) * nx;   # lower row 2 3rd node
          if idone[nod[8]] != 1
            tt0 = sqr2 * slw[iss] + tt[iss];
            tt,ipath,ndid,nmele,idele  = update(tt0,nod[8],tt,ipath,ndid,iss,tmin,
            nmele,idele);
          end
        end

      end

    # order 3, add another 8 nodes
    if j >= 3

      if i >= 2
        nod[9] = i - 1 + (j - 3) * nx;     # upper row 1 1st node
        if idone[nod[9]] != 1
          tt0 = sqr5 * (slw[nod[9]] + slw[nod[1]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[9],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end
      if i <= nx-1
        nod[10] = i + 1 + (j -3) * nx;    # upper row 2nd node
        if idone[nod[10]] != 1
          tt0 = sqr5 * (slw[nod[2]-nx] + slw[nod[2]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[10],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

    end

    if j >= 2

      if i >= 3
        nod[11] = i - 2 + (j - 2) * nx;
        if idone[nod[11]] != 1
          tt0 = sqr5 * (slw[nod[11]] + slw[nod[1]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[11],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end
      if i <= nx-2
        nod[12] = i + 2 + (j - 2) * nx;
        if idone[nod[12]] != 1
          tt0 = sqr5 * (slw[nod[2]] + slw[nod[3]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[12],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

    end

    if j <= nz-1

      if i >= 3
        nod[13] = i - 2 + (j + 0) * nx;
        if idone[nod[13]] != 1
          tt0 = sqr5 * (slw[nod[4]] + slw[nod[4]-1]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[13],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end
      if i <= nx-2
        nod[14] = i + 2 + (j + 0) * nx;
        if idone[nod[14]] != 1
          tt0 = sqr5 * (slw[iss] + slw[nod[5]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[14],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

    end

    if j <= nz-2

      if i >= 2
        nod[15] = i - 1 + (j + 1) * nx;
        if idone[nod[15]] != 1
          tt0 = sqr5 * (slw[nod[4]] + slw[nod[6]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[15],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end
      if i <= nx-1
        nod[16] = i + 1 + (j + 1) * nx;
        if idone[nod[16]] != 1
          tt0 = sqr5 * (slw[iss] + slw[nod[7]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[16],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

    end


    # order 4, add another 16 nodes
    if iorder >= 4

      if i >= 3 && j >= 4
        nod[17] = i - 2 + (j - 4) * nx;
        if idone[nod[17]] != 1
          tt0 = sqr13 * (slw[nod[17]] + 0.5 * slw[nod[17]+nx] + 0.5 * slw[nod[9]]
          + slw[nod[1]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[17],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 2 && j >= 4
        nod[18] = i - 1 + (j - 4) * nx;
        if idone[nod[18]] != 1
          tt0 = sqr10 * (slw[nod[18]] + slw[nod[9]] + slw[nod[1]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[18],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-1 && j >= 4
        nod[19] = i + 1 + (j - 4) * nx;
        if idone[nod[19]] != 1
          tt0 = sqr10 * (slw[nod[18]+1] + slw[nod[9]+1] + slw[nod[2]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[19],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-2 && j >= 4
        nod[20] = i + 2 + (j-4) * nx;
        if idone[nod[20]] != 1
          tt0 = sqr13 * (slw[nod[19]] + 0.5 * slw[nod[10]] + 0.5 * slw[nod[10]-1]
          + slw[nod[2]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[20],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 4 && j >= 3
        nod[21] = i - 3 + (j - 3) * nx;
        if idone[nod[21]] != 1
          tt0 = sqr13 * (slw[nod[21]] + 0.5 * slw[nod[21]+1] + 0.5 * slw[nod[11]]
          + slw[nod[1]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[21],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-3 && j >= 3
        nod[22] = i + 3 + (j - 3) * nx;
        if idone[nod[22]] != 1
          tt0 = sqr13 * (slw[nod[22]-1] + 0.5 * slw[nod[10]] + 0.5 * slw[nod[3]]
          + slw[nod[2]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[22],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 4 && j >= 2
        nod[23] = i - 3 + (j - 2) * nx;
        if idone[nod[23]] != 1
          tt0 = sqr10 * (slw[nod[23]] + slw[nod[11]] + slw[nod[1]]) +tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[23],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-3 && j >= 2
        nod[24] = i + 3 + (j - 2) * nx;
        if idone[nod[24]] != 1
          tt0 = sqr10 * (slw[nod[12]] + slw[nod[3]] + slw[nod[2]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[24],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 4 && j <= nz-1
        nod[25] = i - 3 + (j + 0) * nx;
        if idone[nod[25]] != 1
          tt0 = sqr10 * (slw[nod[25]-nx] + slw[nod[13]-nx] + slw[nod[4]]) +
          tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[25],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-3 && j <= nz-1
        nod[26] = i + 3 + (j + 0) * nx;
        if idone[nod[26]] != 1
          tt0 = sqr10 * (slw[iss] + slw[nod[5]] + slw[nod[5]+1]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[26],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 4 && j <= nz-2
        nod[27] = i - 3 + (j + 1) * nx;
        if idone[nod[27]] != 1
          tt0 = sqr13 * (slw[nod[25]] + 0.5 * slw[nod[13]] + 0.5 *
          slw[nod[13]-nx] + slw[nod[4]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[27],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-3 && j <= nz-2
        nod[28] = i + 3 + (j + 1) * nx;
        if idone[nod[28]] != 1
          tt0 = sqr13 * (slw[iss] + 0.5 * slw[nod[5]] + 0.5 * slw[nod[8]] +
          slw[nod[14]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[28],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i >= 3 && j <= nz-3
        nod[29] = i - 2 + (j + 2) * nx
        if idone[nod[29]] != 1
          tt0 = sqr13 * (slw[nod[27]+1] + 0.5 * slw[nod[13]] + 0.5 * slw[nod[6]]
          + slw[nod[4]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[29],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end


      if i >= 2 && j <= nz-3
        nod[30] = i - 1 + (j + 2) * nx;
        if idone[nod[30]] != 1
          tt0 = sqr10 * (slw[nod[4]] +slw[nod[6]] + slw[nod[15]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[30],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-1 && j <= nz-3
        nod[31] = i + 1 + (j + 2) * nx;
        if idone[nod[31]] != 1
          tt0 = sqr10 * (slw[iss] +slw[nod[7]] +slw[nod[16]-1]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[31],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

      if i <= nx-2 && j <= nz-3
        nod[32] = i + 2 + (j + 2) * nx
        if idone[nod[32]] != 1
          tt0 = sqr13 * (slw[iss] + 0.5 * slw[nod[7]] + 0.5 * slw[nod[8]] +
          slw[nod[16]]) + tt[iss];
          tt,ipath,ndid,nmele,idele  = update(tt0,nod[32],tt,ipath,ndid,iss,tmin,
          nmele,idele);
        end
      end

    end     # order >=4























    # find the next 'source'
    @label loop1
    if nmele[i1st] == 0

      for m in i1st+1: imax
        if nmele[m] != 0
          i1st = m;
          @goto loop2
        end
      end

      @goto sorend

    @label loop2
      for mm in 1:nmele[i1st]
        idone[idele[mm,i1st]] = 1;   # all selected to be sources
      end

    end

    @label loop3
    iss = idele[nmele[i1st],i1st];
    nmele[i1st] = nmele[i1st] - 1;
    ndid = ndid - 1;

    if indcn[iss] == 1
      if nmele[i1st] > 0
        @goto loop3
      end
      if nmele[i1st] == 0
        @goto loop1
      end
    end

    if ndid == 0
      @goto sorend
    end
    i,j = denum(iss,nx,nz);

    @goto order     # start the second time calculation

    @label sorend

    # extract raypaths，if no raypaths request, just skip this paragraph
    for ir in 1: nr[is]

      ix = Int64(floor((rx[ir,is] - x0) / dx + 1));
      iz = Int64(floor((rz[ir,is] - z0) / dx + 1));
      iirc = Int64(floor(ix + (iz - 1) * nx));
      icnt = 0;
    @label getsor
      iirp = ipath[iirc];
      if iirp == 0
        @goto etrend
      end
      inod,nd = interp0(iirc,iirp,nd,dx,nx,nz);  # ??why array should be in, Int cannot
      for k in 1:inod
        icnt = icnt + 1;
        iasen[ip+icnt] = nd[k];
      end

      iirc = iirp;
      @goto getsor
    @label etrend
      iasen[ip] = icnt;
      ip = ip + icnt + 1;
    end

    if ip > msmax*ns
      msg = "ip > mamax in raytracing";  # add a return here
    end

    # get the traveltime of each receiver
    ttime = timrec(is,nr,nx,x0,z0,dx,rx,rz,tt,ttime);

  end     # !! source loop end

  return tt,ttime
end
