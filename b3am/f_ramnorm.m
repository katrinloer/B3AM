function data_norm = f_ramnorm(data,fmin,sr)
% Note: relative amplitudes of components are not preserved - not suitable
% for ellipticity analysis!

ldata = length(data);

% Running-absolute-mean normalisation after Benson et al. (2007)
Tmax = 1/fmin;
lwin = round(0.5*Tmax*sr);
if mod(lwin,2)==0
    lwin = lwin+1;
end
N = (lwin-1)/2;

data_help = zeros(ldata+2*N,1);
data_help(N+1:end-N) = data;

data_norm = zeros(size(data));
c = 0;
for i = (N+1):(N+ldata)
    c = c+1;
    w = mean(abs(data_help(i-N:i+N)));
    data_norm(c) = data_help(i)./w;
end

end