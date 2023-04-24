%% Write station_info.txt file for B3AM from DAT

% BEFORE YOU RUN THIS SCRIPT:
% - load DAT file that contains DAT.h.coords into workspace
% - define name for new station file (incl path) in l. 7

statname = './station_utm_BH_2002093.txt';
fid = fopen(statname,'w');
nstat = length(DAT(1).h.coords);
for i = 1:nstat
    sname = DAT.h.stations{i};
    x = DAT.h.coords(i,1);
    y = DAT.h.coords(i,2);
    fprintf(fid,'%s %f %f\n',sname,x,y);
end
fclose(fid); 