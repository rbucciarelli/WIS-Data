function [sp1d] = sp2d_to_1d(sp2d,dtheta)
f=fscanf(fid,'%g',NF);              %- freqs
dir=fscanf(fid,'%g',ND);            %- direction in radians
fprop=f(2)/f(1);                    %- freq fraction
fdum=[f(1)/fprop;f];                %- place holder/normalization
df=diff(fdum);                      %- bandwidth
dtheta=abs(dir(2)-dir(1));          %- direction bin size (radians)


end
