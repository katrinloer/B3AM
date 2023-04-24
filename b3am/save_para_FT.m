function DFT = save_para_FT(DATname,data,procpars,t0,coords,stations,dt_new,f,f0)

% Extract spectral data for a given frequency
DFT.data = data;

% Store Fourier processing parameters in a header
DFT.procpars = procpars;

% Add frequency of this spectral 'layer' to header
DFT.h.f0 = f0;

% Add header information (start time, dt, coordinates, etc.)
DFT.h.t0 = t0;
DFT.h.coords = coords;
DFT.h.stations = stations;
DFT.h.dt = dt_new;
DFT.h.f = f; % Also add the full frequency vector

% Save
save(DATname,'DFT')%,'-v7.3');

end