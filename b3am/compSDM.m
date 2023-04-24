function [fnew,S] = compSDM(f,DFTS3C,NBlock,doNorm)
% !!! TEST: CHECK LINE 120 !!!
% Compute SDM of a batch of spectrograms (f,t) and perform block-averaging
% as well as frequency stacking (reducing frequency sampling in the
% process, i.e. no smoothing).
% 
% If DFTS3C is a gpuArray object, the output S,Sn should also be a gpuArray
% object.
% 
% 
% USAGE ======================
% 
%   [fnew,S] = compSDM(f,DFTS3C,NBlock,doNorm)
%
% 
% INPUT ======================
% f             (F,1) real. Frequency vector
% DFTS3C        (F,L,3*n) complex. Spectrogram (F freq. bins, L time windows)
%               for all 3*n receivers in the array.
% NBlock        Int. Number of time windows to average SDM over
% doNorm        Logical. Normalize Fourier amplitudes to 1 before computing
%               SDM. Note that in the current implementation this will
%               destroy the ellipticity information for Rayleigh waves.
%               
% 
% OUTPUT ======================
% fnew          (F) real. New frequency vector (frequency averaging is
%               deactivated as of now. fnew is identical to f.
% S             (3*n,3*n,F,L) complex. Spectral density matrix between the
%               3*n receivers at each of the F frequency bins (may be just
%               1) and L time windows.
% 
% 
% 
% NOTES ======================
% This function is part of a package for three-component wavenumber
% analysis. A possible workflow to use this package comes in the form of
% two sequential scripts that can be run.
% 
% 
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
% 
% Early version, please don't distribute yet.
% 


% Set normalization off by default
if nargin<5
    doNorm = false;
end


% Frequency averaging, must be a natural number. It defines the new
% frequency resolution where each nfstk*df Hz there will be an average of
% nfstk original frequency bins
nfstk = 1; 

% Block-averaging over time, again a natural number
ntstk = NBlock;

if nfstk>1
    nfreq = floor(length(f)/nfstk);
    inds = reshape(1:nfreq*nfstk,nfstk,nfreq);
    
    % Construct new frequency vector
    fnew = mean(f(inds),1)';
else
    fnew = f;
end

% Get M(number of freqs), N(number of time segs), K(number of locations)
% Spectrogram (MxN) synchronicity over all locations is assumed
[M,N,K] = size(DFTS3C);

% Each column is the full spectrogram (freq x time) of a sensor
% location-component combination
DFTS3C = reshape(DFTS3C, M*N , K );

% Smoothing kernel for block-averaging
smkernel = ones(1,ntstk);
smkernel = smkernel/sum(smkernel(:));


%% Compute SDM among sensors

% fprintf('Start SDM computation.\n');

S = single(zeros(K,K,length(fnew),N));

% Compute lower triangle SDM (and put conj(.) to upper simultaneously
for j = 1:K
    for i = j:K
        
        %% Compute SDM
        % Compute (i,j) entry for all freq/times
        if doNorm
            tmp = (DFTS3C(:,j)./abs(DFTS3C(:,j))) .* (conj(DFTS3C(:,i))./abs(DFTS3C(:,i)));
        else
            tmp = (DFTS3C(:,j)) .* (conj(DFTS3C(:,i)));
        end
        
        % Reshape into spectrogram of (i,j) pair
        tmp = reshape(tmp,M,N);
        
        if nfstk>1
            % Stack nfstk frequencies and write back to same matrix (works
            % because small indices are not used anymore)
            for k = 1:size(inds,2)
                tmp(k,:) = mean(tmp(inds(:,k),:),1);
            end
            % Remove upper indices which are not needed anymore
            tmp = tmp(1:nfreq,:);
        end
        
        % Smooth in time (= block averaging)
        tmp = smoothmat( tmp,smkernel );
        
        % NEW!!!
        tmp = conj(tmp); % KL: This is necessary as the above compuation swaps i and j
        
        % Store that spectrogram and its compl. conj, at (i,j) and (j,i)
        S(i,j,:,:) = tmp;
        S(j,i,:,:) = conj(tmp);
    
    
    end
end

% Clear Fourier amplitudes, they are no longer needed
clear DFTS3C

% fprintf('SDM computation completed.\n');

end


function mat = smoothmat(mat,smkernel)
% 
% Smooth a matrix along its rows and/or columns. Essentially a conv2
% operation, but with appropriate cropping to keep the size of the data
% matrix the same. Smoothing window (1xn, nx1, nxm) must be supplied.
% 
%   mat = nri_smoothmat(mat,smkernel)
% 
% Typically:
% smkernel = ones(1,nsmo)/nsmo;
% OR
% hwin = hanning(nsmo);
% smkernel = hwin * hwin';
% smkernel = smkernel/sum(smkernel(:));
% OR similar...
% 


[nsmoC,nsmoR] = size(smkernel);
if min(nsmoC,nsmoR)~=0
    if nsmoC==1 % row smooth
        mat = conv2(mat,smkernel);  % Order of arguments of conv2 does not matter
        mat = mat(:,1+(nsmoR-1)/2:end-(nsmoR-1)/2);
    elseif nsmoR==1
        mat = conv2(mat,smkernel);  
        mat = mat(1+(nsmoC-1)/2:end-(nsmoC-1)/2,:);
    else
        mat = conv2(mat,smkernel);  
        mat = mat(1+(nsmoC-1)/2:end-(nsmoC-1)/2 , 1+(nsmoR-1)/2:end-(nsmoR-1)/2);
    end
end



end
