function save_para_amp(sname,kr_all,kth_all,pola_all,pola_ind,a_all,wave_ind)

save(sname,'kr_all','kth_all','pola_all','pola_ind','a_all','wave_ind');
% kr_all: horizontal wavenumber
% kth_all: azimuth
% pola_all: polarization parameters (4)
% pola_ind: polarization index
% a_all: amplitude of max. beam response
% wave_ind: wavetype index (0 P, 1=SV, 2=SH/Love, 3=retr. Rayleigh, 4=pro.
% Rayleigh)

end