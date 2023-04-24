function data_norm = f_ramnorm3C(data,fmin,sr)

ldata = length(data);

% Running-absolute-mean normalisation after Benson et al. (2007)
Tmax = 1/fmin;
lwin = round(0.5*Tmax*sr);
if mod(lwin,2)==0
    lwin = lwin+1;
end
N = (lwin-1)/2;

data_E = zeros(ldata+2*N,1);
data_E(N+1:end-N) = data(:,1);

data_N = zeros(ldata+2*N,1);
data_N(N+1:end-N) = data(:,2);

data_Z = zeros(ldata+2*N,1);
data_Z(N+1:end-N) = data(:,3);

data_norm = zeros(size(data));
c = 0;
for i = (N+1):(N+ldata)
    
    c = c+1;
    
    % Normalize with respect to z-component
    wZ = mean(abs(data_Z(i-N:i+N)));
    
    data_norm(c,1) = data_E(i)./wZ;
    data_norm(c,2) = data_N(i)./wZ;
    data_norm(c,3) = data_Z(i)./wZ;
end

end