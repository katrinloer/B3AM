function [P, Q, kr, kth, polstates] = f_FK3C_fast(DFT,procpars)

DFT3C = DFT.data;

% Get M(number of freqs), N(number of time segs), K(number of locations)
% Spectrogram (MxN) synchronicity over all locations is assumed
[Nf,Nt,K] = size(DFT3C);

% Each column is the full spectrogram (freq x time) of a sensor
% location-component combination
DFT3C = reshape(DFT3C, Nf*Nt , K );

% Get the dimensions (K: # receivers, Nf: # of frequencies, Nt: # of time
% windows).
% [~,K,Nf,Nt] = size(S); %fprintf('# of time windows: %d',Nt)
if Nf>1
    error('This version of compFK3C currently only supports single frequency computation.');
end

%% Block-averaging over time

ntstk = procpars.Nblock;

% Smoothing kernel for block-averaging
smkernel = ones(1,ntstk);
smkernel = smkernel/sum(smkernel(:));

for i = 1:K
    
    % Smooth in time (= block averaging)
    DFT3C(:,i) = smoothmat( DFT3C(:,i),smkernel );

%     tmp = DFT3C(:,i);
% 
%     % Smooth in time (= block averaging)
%     tmp = smoothmat( tmp,smkernel );
%     
%     sDFT(:,i) = tmp;

end

% Subsample in time
%-------------------
% It is reasonable to subsample at about half the block-averaging
% duration
ind = 1:procpars.ntsubsample:Nt;
DFT3C = DFT3C(ind,:);

%% Compute 3C wavenumber spectra
%-------------------------------

[P,Q,kr,kth,polstates] = compFK3C22_fast(DFT.h.coords,DFT3C,...
    procpars.kgrid,...
    procpars.agrid,...
    procpars.polstates,...
    procpars.method);
% [P,Q,kr,kth,EVals,polstates] = compFK3C22_fast(DFT.h.coords,DFT.data,...
%     procpars.kgrid,...
%     procpars.method,...
%     procpars.dopolar);
end
