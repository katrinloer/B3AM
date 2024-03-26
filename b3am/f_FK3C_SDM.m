function [P, Q, kr, kth, polstates] = f_FK3C_SDM(DFT,procpars)

% Construct time vector (start time plus signal duration converted to
% MATLAB days)
t0 = DFT.h.t0(1);
t = t0 + (0:size(DFT.data,2)-1)*DFT.h.dt /24/3600;

% Compute a series of SDMs for each time window, including
% block-averaging (effectively a moving average)
% S is a (K*3, K*3, 1, NanalysisWindows) structure. S(:,:,i) is the spectral density matrix of the i-th analysis period.
[f0,S] = compSDM(DFT.h.f0,...
    reshape(DFT.data,1,size(DFT.data,2),size(DFT.data,3)),...
    procpars.Nblock);

% Subsample in time
%-------------------
% It is reasonable to subsample at about half the block-averaging
% duration
ind = 1:procpars.ntsubsample:length(t);
t = t(ind);
S = S(:,:,:,ind);

% Compute 3C wavenumber spectra
%-------------------------------
[P,Q,kr,kth,EVals,polstates] = compFK3C22(DFT.h.coords,S,...
    procpars.kgrid,...
    procpars.agrid,...
    procpars.polstates,...
    procpars.method);

end
