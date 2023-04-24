function data_one = f_onebit(data)

% One-bit normalization
ldata = length(data);
data_one = zeros(ldata,1);
for l = 1:ldata
    if data(l) > 0
        data_one(l) = 1;
    else
        data_one(l) = -1;
    end
end

end