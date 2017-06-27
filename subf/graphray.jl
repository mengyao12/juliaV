function graphray(nx::Int64,nz::Int64,nm::Int64,slw::Array{Float64,1},
  izdn::Int64,itop::Array{Int64,1},ixr1::Array{Int64,1},
  ixr2::Array{Int64,1},ns::Int64,dx::Float64,imax::Int64,
  x0::Float64,z0::Float64,sx::Array{Float64,1},sz::Array{Float64,1},
  nxz::Int64)

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
  i1st = Int64;
  ndid = Int64;
  i = Int64;
  j = Int64;
  iss = Int64;


  for ii in 1:nm
    idone[ii] = 1;
  end

  for is in 1:10   # !! source loop

    # initialize the grid for each source
    idone,tt = washtomo(is,ns,nx,nm,ixr1,ixr2,itop,ibt,idone,tt);
    nmele = washnm(nmele,imax);
    indcn = washnm(indcn,nm);

    # calculate the tt of source
    nmele,idele,ipath,idone,tt,i1st,ndid,iss,i,j = timsor(nx,x0,z0,dx,is,sx,sz,
    nmele,idele,nm,ipath,slw,idone,tt);

    # order 2, 8 nodes
    



  end

  return nmele,idele,ipath,idone,tt,i1st,ndid,iss,i,j
end
