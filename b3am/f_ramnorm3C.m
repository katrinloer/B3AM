function [dataE_norm,dataN_norm,dataZ_norm] = f_ramnorm3C(dataE,dataN,dataZ,fmin,sr)
% Preserves relative amplitudes of components in each time window
% Note: as a result, if a spike occurs on one component only, this will not be removed

ldata = length(dataE);

% Running-absolute-mean normalisation after Benson et al. (2007)
Tmax = 1/fmin;
lwin = round(0.5*Tmax*sr);
if mod(lwin,2)==0
    lwin = lwin+1;
end
N = (lwin-1)/2;

data_E = zeros(ldata+2*N,1);
data_E(N+1:end-N) = dataE;

data_N = zeros(ldata+2*N,1);
data_N(N+1:end-N) = dataN;

data_Z = zeros(ldata+2*N,1);
data_Z(N+1:end-N) = dataZ;

dataE_norm = zeros(size(dataE));
dataN_norm = zeros(size(dataN));
dataZ_norm = zeros(size(dataZ));
c = 0;
for i = (N+1):(N+ldata)
    
    c = c+1;
    
    % Normalize with respect to z-component
    wZ = mean(abs(data_Z(i-N:i+N)));
    
    dataE_norm(c) = data_E(i)./wZ;
    dataN_norm(c) = data_N(i)./wZ;
    dataZ_norm(c) = data_Z(i)./wZ;
end

end