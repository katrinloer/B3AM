%% Plot beam response P and polarisation matrix Q for single frequency and time window
% Note: variables need to be in the workspace (e.g., after running b3am.m)

% close all

savefig = 0;    % save figures? 0 (no) or 1 (yes)
nonoise = 0;    % plot beam levels above noise (3x standard deviation)
Qmode = 'wave'; % 'pola' for polarisation indices, 'wave' for wave types
fsize = 14;     % Figure fontsize

addpath('/Users/kloer/Documents/MATLAB/crameri');
addpath('/Users/kloer/Documents/PROJECTS/FKanalysis/B3AM/plot');

figure('Color','white','Position',[50 150 297*3 210*3]);

for itime = 8%1:1:size(P,3) % choose time window(s)

ix = 1 ; % largest maximum = 1

% New
wave_ix = wave_ind(itime,ix); 

% Wave type
if wave_ix == 0
    w0 = 'P wave';
elseif wave_ix == 1
    w0 = 'SV wave';
elseif wave_ix == 2
    w0 = 'Love/SH wave';
elseif wave_ix == 3
    w0 = 'Retrograde Rayleigh wave';
elseif wave_ix == 4
    w0 = 'Prograde Rayleigh wave';
end

% Azimuth 
azi = kth_max(itime,ix)*180/pi;

% Dip
dip = pola_max(itime,ix,2) * 180/pi;

% Ellipticity
ell = pola_max(itime,ix,3);

% Velocity
vel = f/kr_max(itime,ix) * sin(pola_max(itime,ix,2));

%% Plot
    
tl = tiledlayout(1,2);
title(tl, sprintf('%s: azimuth = %3.0f°, dip = %3.0f°, ellipticity = %2.1f', w0, azi, dip, ell), 'Fontsize', 18);

Pplot = P(:,:,itime);

% Beam response
ax1 = nexttile;
imagesc(kr*1000,kth*180/pi,Pplot)
xlabel('wavenumber in 1/km')
ylabel('azimuth in °')
title('Beam response')
cb = colorbar;
cb.Label.String = 'beampower';
cb.Label.FontSize = fsize;
cmap = crameri('batlow');
colormap(ax1, cmap);
if nonoise
    Pmax = max(Pplot(:));%max(max(Pplot));
    Pmean = mean(Pplot(:));%mean(mean(Pplot));
    Pstd = std(Pplot(:));%std(std(Pplot));
    caxis([Pmean+3*Pstd Pmax])
end
set(gca,'Fontsize',fsize)

% Plot marker for maximum
hold on
plot3(kr_max(itime,ix)*1000,kth_max(itime,ix)*180/pi,max(max(Pplot))+1, 'w+','MarkerSize',12)

% Polarisation...
Qplot = Q(:,:,itime);
ax2 = nexttile;

% ... use polarisation index (rows in polstates)
if strcmp(Qmode,'pola')

    imagesc(kr*1000,kth*180/pi,Qplot)
    xlabel('wavenumber in 1/km')
    ylabel('azimuth in °')
    title('Polarisation')
    cb = colorbar;
    cb.Label.String = 'polariation ID';
    cb.Label.FontSize = fsize;
    cmap = crameri('batlow');
    colormap(ax2,cmap);
    set(gca,'Fontsize',fsize)
    hold on

% ...use wave type index (0-4)
elseif strcmp(Qmode,'wave')

    Qw = zeros(size(Qplot));

    for i = 1:size(Qplot,1)
        for j = 1:size(Qplot,2)

            pola_help = polstates(Qplot(i,j),:);

            if pola_help(3) == 0 && pola_help(4) == pi
                % P-wave
                Qw(i,j) = 0;
            elseif pola_help(3) == 2 && pola_help(4) == pi
                % SV-wave
                Qw(i,j) = 1;
            elseif pola_help(3) == 2 && pola_help(4) == pi/2
                % Love wave
                Qw(i,j) = 2;
            elseif pola_help(3) > 0 && pola_help(3) < 2 && pola_help(4) == 0
                % Retrograde Rayleigh wave
                Qw(i,j) = 3;
            elseif pola_help(3) > 0 && pola_help(3) < 2 && pola_help(4) == pi
                % Prograde Rayleigh wave
                Qw(i,j) = 4;
            end

        end
    end

    imQ = imagesc(kr*1000,kth*180/pi,Qw);
    xlabel('wavenumber in 1/km')
    ylabel('azimuth in °')
    title('Wave type')
    cb = colorbar;
    cb.Ticks = 0:1:4;
    cb.TickLabels = {'P','SV','SH','rR','pR'};
    load batlowS.mat
    cmapQ = batlowS(1:5,:);
    colormap(ax2, cmapQ);
    caxis([0 4])
    set(gca,'Fontsize',fsize)
    hold on

end

% Plot marker for maximum
hold on
plot3(kr_max(itime,ix)*1000,kth_max(itime,ix)*180/pi,max(max(Pplot))+1, 'w+','MarkerSize',12)

%% Save figure

if savefig
    figname = sprintf('%sFigures/b3am_twin%03d.jpg',indir,itime);
    print(figname,'-djpeg')
end

pause(1)

end

%% EOF