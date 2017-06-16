function readparm(parmpath::AbstractString)
  f = open(parmpath);
  s = readlines(f);
  ttpth = s[2];
  ttpth = ttpth[1:end-1];
  mdlpth = s[4];
  mdlpth = mdlpth[1:end-1];
  outpth = s[6];
  outpth = outpth[1:end-1];
  maxitr = parse(Int64,s[8]);
  maxV = parse(Float64,s[10]);
  minV = parse(Float64,s[12]);
  tau = parse(Float64,s[14]);
  hvs = parse(Float64,s[16]);
  close(f)

  return ttpth,mdlpth,outpth,maxitr,maxV,minV,tau,hvs

end
