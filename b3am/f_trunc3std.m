function data = f_trunc3std(data)

% Define threshold as 3x standard deviation
maxdata = mean(data) + 3*std(data);
mindata = mean(data) - 3*std(data);

% Truncate data at threshold
data(data > maxdata) = maxdata;
data(data < mindata) = mindata;

end