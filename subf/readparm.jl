function readparm(parmpath::AbstractString)
f = open(parmpath);
s = readlines(f);
ttpth = s[2];
mdlpth = s[4];
outpth = s[6];
maxitr = s[8];
maxV = s[10];
minV = s[12];
tau = s[14];
hvs = s[16];
return ttpth,mdlpth,outpth,maxitr,maxV,minV,tau,hvs
close(f)
end
