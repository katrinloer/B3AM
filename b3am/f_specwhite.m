%% Spectral whitening

function data_white = f_specwhite(data,fs)

L = length(data);
n = L;%2^nextpow2(L); % back transformation not working if n ~= L

% Fourier transform (time -> freq)
D = fft(data,n);

f = fs*(0:(n/2))/n;
D_abs = abs(D/n);

% Envelope
nav = 10;%100; % number of elements to average over
[env, ~] = envelope(D_abs,nav,'rms');

% Divide by envelope
D_white = D ./ env;

% Fourier transform (freq -> time)
data_white = real(ifft(D_white));

end