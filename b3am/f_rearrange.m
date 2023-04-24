function data_re = f_rearrange(data,nstat,rstep)

for i = 1:rstep
    
    data_re(:,nstat*(i-1)+1:nstat*i) = data(:,i:rstep:end);
    
end

end

% Resorts data/station/channel at location
% 1,4,7,... to 1,2,3,...,n          (channel E)
% 2,5,8,... to n+1,n+2,n+3,...2n    (channel N)
% 3,6,9,... to 2n+1,2n+2,2n+3,...3n (channel Z)