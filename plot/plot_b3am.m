%% Plot functions for b3am
% wavefield composition
% surface wave f-k histograms and dispersion curves
% direction of arrival
%
% Katrin Loer
% Nov 2016
% Last modified Mar 2022: implementing accumarray and findpeak
%--------------------------------------------------------------------------

clear
close all

% Add paths to figure printing functions and colourmaps
addpath('/Users/s06kl9/Projects/FKanalysis/SCRIPTS/export_fig');
addpath('/Users/s06kl9/Projects/FKanalysis/SCRIPTS/crameri_v1/crameri/');
addpath('/Users/s06kl9/Projects/FKanalysis/SCRIPTS/crameri_v1/CategoricalPalettes/');

% INPUT: Set path to B3AM output (kmax files)
dir_in =  '/Users/s06kl9/Projects/FKanalysis//GEMEX/DATA/2017/LH/DAT_B3AM_TEST/kmax/';

% OUTPUT: Define path for output figures
dir_out = '/Users/s06kl9/Projects/FKanalysis/B3AM/Figures';

% Choose plot options:
plotcomp = true;           % plot wavefield composition in bar plot
plotcompcurve = true;      % plot wavefield composition as curves
plotkhisto = true;         % plot wavenumber histogram
plotvdisp = true;          % plot dispersion curves as velocity vs frequency
plotDOA = true;

% Save as figure?
savecomp = false;
savecompcurve = false;
savekhisto = false;
savevdisp = false;
saveDOA = false;

% Minimum signal to noise ratio in histograms
SNR = 20; % amplitude/mean(amplitudes) for each frequency --> threshold for maxima in histogram; 20 works well for Groningen ('MAX1')

maxflag = 'MAX1';%'NOMAX' or 'MAX1'

countflag = 'amp'; % 'amp' or 'noamp' (consider amplitudes when counting detections or not - for wavefield composition plot only)

%% Don't change:

% Load processing parameters
load([dir_in 'procpars.mat'])

freqs = procpars.freqs;
kmin = min(procpars.kgrid);
kmax = max(procpars.kgrid);

% Load categorical colour palette
load batlowS.mat % for line plotting

infiles = dir([dir_in,'kmax_*.mat']);    
indate = infiles(1).name(9:15);
nw = infiles(1).name(6:7);

nfreq = length(freqs); % number of frequencies (see infiles)
ndays_total = length(infiles)/nfreq;
ndays_stack = ndays_total;
daystep = 1;

% Allocate cell arrays
kr_retro_all = cell(length(freqs),1);
kr_pro_all = cell(length(freqs),1);
kr_love_all = cell(length(freqs),1);
kr_p_all = cell(length(freqs),1);
kr_sv_all = cell(length(freqs),1);

kth_retro_all = cell(length(freqs),1);
kth_pro_all = cell(length(freqs),1);
kth_love_all = cell(length(freqs),1);
kth_p_all = cell(length(freqs),1);
kth_sv_all = cell(length(freqs),1);

a_retro_all = cell(length(freqs),1);
a_pro_all = cell(length(freqs),1);
a_love_all = cell(length(freqs),1);
a_p_all = cell(length(freqs),1);
a_sv_all = cell(length(freqs),1);

e_retro_all = cell(length(freqs),1);
e_pro_all = cell(length(freqs),1);
dip_p_all = cell(length(freqs),1);
dip_sv_all = cell(length(freqs),1);

%% WAVE TYPE ANALYSIS

% Loop over frequencies
for ff = 1:nfreq

    f = freqs(ff);%str2double(infiles(ff).name(end-4:end));
    
    % Loop over total number of days
    for tt = 1:daystep:ndays_total
        
        if ndays_total-tt+1 < ndays_stack
            continue
        end
        
        % Count surface and body waves in frequency bin
        
        c_pro = 0; kr_pro = 0; kth_pro = 0; w_pro = 0;
        c_retro = 0; kr_retro = 0; kth_retro = 0; w_retro = 0;
        c_love = 0; kr_love = 0; kth_love = 0; w_love = 0;
        c_p = 0; kr_p = 0; kth_p = 0; w_p = 0;
        c_sv = 0; kr_sv = 0; kth_sv = 0; w_sv = 0;
        
        % Loop over days to stack
        for dd = tt:min(tt+ndays_stack-1,ndays_total)
            
            datafile = [dir_in, infiles(ff+(dd-1)*nfreq).name];
            if dd == tt
                outfile_help = datafile;
            end
            load(datafile)
            
            if strcmp(maxflag,'MAX1')
                nmax = 1;
            elseif strcmp(maxflag,'MAX3')
                nmax = 3;
            elseif strcmp(maxflag,'MAX5')
                nmax = 5;
            elseif strcmp(maxflag,'NOMAX')
                nmax = size(kr_all,2);
            end
            
            L = size(kr_all,1);
            
            for n = 1:nmax % number of maxima to count
                for j = 1:L % time windows
                    
                    if kr_all(j,n) <= kmin || kr_all(j,n) >= kmax
                        continue
                    end
                    
                    % P-wave
                    if wave_ind(j,n) == 0
%                     if pola_ind(j,n) > 1 && pola_ind(j,n) <= 10
                        c_p = c_p + 1;
                        kr_p(c_p) = kr_all(j,n);
                        kth_p(c_p) = kth_all(j,n);
                        w_p(c_p) = a_all(j,n); % amplitude/weighting factor
                        dip_p(c_p) = pola_all(j,n,2);
                    % Love wave
                    elseif wave_ind(j,n) == 2
%                     elseif pola_ind(j,n) == 11
                        c_love = c_love + 1;
                        kr_love(c_love) = kr_all(j,n);
                        kth_love(c_love) = kth_all(j,n);
                        w_love(c_love) = a_all(j,n);
                    % SV-wave
                    elseif wave_ind(j,n) == 1
%                     elseif pola_ind(j,n) > 11 && pola_ind(j,n) <= 21
                        c_sv = c_sv + 1;
                        kr_sv(c_sv) = kr_all(j,n);
                        kth_sv(c_sv) = kth_all(j,n);
                        w_sv(c_sv) = a_all(j,n);
                        dip_sv(c_sv) = pola_all(j,n,2);
                    % Prograde Rayleigh wave
                    elseif wave_ind(j,n) == 4
%                     elseif pola_ind(j,n) > 21 && pola_ind(j,n) <= 40
                        c_pro = c_pro + 1;
                        kr_pro(c_pro) = kr_all(j,n);
                        kth_pro(c_pro) = kth_all(j,n);
                        w_pro(c_pro) = a_all(j,n);
                        e_pro(c_pro) = pola_all(j,n,3);
                    % Retrograde Rayleigh wave
                    elseif wave_ind(j,n) == 3
%                     elseif pola_ind(j,n) > 40 && pola_ind(j,n) <= 59
                        c_retro = c_retro + 1;
                        kr_retro(c_retro) = kr_all(j,n);
                        kth_retro(c_retro) = kth_all(j,n);
                        w_retro(c_retro) = a_all(j,n);
                        e_retro(c_retro) = pola_all(j,n,3);
                    end
                    
                end
            end
        end
        
        % Number of wave types
        nwaves(ff,1:5) = [c_retro, c_pro, c_love, c_sv, c_p];
        
        % Amplitudes of wave types
        awaves(ff,1:5) = [sum(w_retro),sum(w_pro),sum(w_love),sum(w_sv),sum(w_p)];
        
        % Wavenumber
        if c_retro > 0; kr_retro_all{ff} = kr_retro; kr_retro_mean(ff) = mean(kr_retro); end
        if c_pro > 0; kr_pro_all{ff} = kr_pro; kr_pro_mean(ff) = mean(kr_pro); end
        if c_love > 0; kr_love_all{ff} = kr_love; kr_love_mean(ff) = mean(kr_love); end
        if c_sv > 0; kr_sv_all{ff} = kr_sv; kr_sv_mean(ff) = mean(kr_sv); end
        if c_p > 0; kr_p_all{ff} = kr_p; kr_p_mean(ff) = mean(kr_p); end
        
        % Ellipticities
        if c_retro > 0; e_retro_all{ff} = e_retro; eretro_m(ff) = mean(e_retro); end
        if c_pro > 0; e_pro_all{ff} = e_pro;  epro_m(ff) = mean(e_pro); end
        
        % Direction of arrival (DOA)
        if c_retro > 0; kth_retro_all{ff} = kth_retro; kth_retro_mean(ff) = mean(kth_retro); end
        if c_pro > 0; kth_pro_all{ff} = kth_pro; kth_pro_mean(ff) = mean(kth_pro); end
        if c_love > 0; kth_love_all{ff} = kth_love; kth_love_mean(ff) = mean(kth_love); end
        if c_sv > 0; kth_sv_all{ff} = kth_sv; kth_sv_mean(ff) = mean(kth_sv); end
        if c_p > 0; kth_p_all{ff} = kth_p; kth_p_mean(ff) = mean(kth_p); end
        
        % Dip for body waves
        if c_sv > 0; dip_sv_all{ff} = dip_sv; dip_sv_mean(ff) = mean(dip_sv); end
        if c_p > 0; dip_p_all{ff} = dip_p; dip_p_mean(ff) = mean(dip_p); end
        
        % Amplitudes
        if c_retro > 0; a_retro_all{ff} = w_retro; end
        if c_pro > 0; a_pro_all{ff} = w_pro;  end
        if c_love > 0; a_love_all{ff} = w_love;  end
        if c_sv > 0; a_sv_all{ff} = w_sv; end
        if c_p > 0; a_p_all{ff} = w_p; end
    end
    
    
end

%% Plot wavefield composition

if plotcomp

pfreqs = freqs;

if strcmp(countflag,'amp')    
    pwaves = awaves;
    ystr_rel = 'relative amplitudes in %';
    ystr_abs = 'absolute amplitudes';
    tstr_rel = 'Wavefield composition: relative (amplitudes)';
    tstr_abs = 'Wavefield composition: absolute (amplitudes)';
elseif strcmp(countflag,'noamp')
    pwaves = nwaves;
    ystr_rel = 'relative counts in %';
    ystr_abs = 'absolute counts';
    tstr_rel = 'Wavefield composition: relative (counts)';
    tstr_abs = 'Wavefield composition: absolute (counts)';
end

% Relative counts
%----------------
figure('Color','white','Position',[2000 10 1000 500]);
% colormap(batlowS);
b = bar(pfreqs,pwaves./sum(pwaves,2).*100,1,'stacked','FaceColor','flat');
for k = 1:5
    b(k).FaceColor = batlowS(k,:);
end
legend('Rayleigh (retro)','Rayleigh (pro)','Love/SH','SV','P','Location','northeast')
set(gca,'xtick',0.2:0.1:1.2,'Fontsize',18)
xlabel('frequency in Hz')
ylabel(ystr_rel)
axis tight 
title(tstr_rel)

if savecomp
    fname = sprintf('%s/comp_%s_rel_%s_%s_%s.png',...
        dir_out,countflag,maxflag,nw,indate);
    export_fig(fname,'-png');
end

% Absolute counts
%----------------
figure('Color','white','Position',[2000 10 1000 500]);
b2 = bar(pfreqs,pwaves,1,'stacked','FaceColor','flat');
for k = 1:5
    b2(k).FaceColor = batlowS(k,:);
end
legend('Rayleigh (retro)','Rayleigh (pro)','Love/SH','SV','P','Location','northeast')
% set(gca,'xtick',freqs(3:4:end),'Fontsize',18)
set(gca,'xtick',0.2:0.1:1.2,'Fontsize',18)
xlabel('frequency in Hz')
ylabel(ystr_abs)
axis tight 
title(tstr_abs)
% ylim([0 0.15*10^5])

if savecomp
    fname = sprintf('%s/comp_%s_abs_%s_%s_%s.png',...
        dir_out,countflag,maxflag,nw,indate);
    export_fig(fname,'-png');
end

% pause(1)

end

%% Plot individual wave contributions (as curves, not bars)

if plotcompcurve

if strcmp(countflag,'noamp')
    
% Counts
figure('Color','white','Position',[2000 10 1000 500]);
title('Wavefield composition: counts')

hold on
p1 = plot(freqs,nwaves,'LineWidth',2);
for k = 1:5
    p1(k).Color = batlowS(k,:);
end

legend('Rayleigh (retro)','Rayleigh (pro)','Love/SH','SV','P','Location','northeast')

box on
grid on
xlabel('frequency in Hz')
ylabel('absolute counts')
set(gca,'Fontsize',18)

if savecompcurve
    fname = sprintf('%s/compcurve_noamp_%s_%s_%s.png',...
        dir_out,maxflag,nw,indate);
    export_fig(fname,'-png');
end

elseif strcmp(countflag,'amp')
    
% Amplitudes
figure('Color','white','Position',[2000 10 1000 500]);
title('Wavefield composition: absolute (amplitudes)')

hold on
p2 = plot(freqs,awaves,'LineWidth',2);
for k = 1:5
    p2(k).Color = batlowS(k,:);
end
legend('re. Rayleigh', 'pr. Rayleigh', 'Love', 'SV', 'P')

box on
grid on
xlabel('frequency in Hz')
ylabel('cummulative amplitude')
set(gca,'Fontsize',18)

if savecompcurve
    fname = sprintf('%s/compcurve_amp_%s_%s_%s.png',...
        dir_out,maxflag,nw,indate);
    export_fig(fname,'-png');
end

end

end

%% COMPUTE WAVENUMBER HISTOGRAMS

% Wavenumber vs frequency
kbins = procpars.kgrid;%linspace(kmin,kmax,201);
kstep = kbins(2)-kbins(1);
% kgrid = linspace(kmin,kmax,201);
% kstep = kgrid(2)-kgrid(1);
% kbins = linspace(kmin-kstep/2,kmax+kstep/2,202);

for i = 1:size(freqs,2)
    
    % Wavenumber histograms
    [Nretro,binretro,idxretro]  = histcounts(kr_retro_all{i},kbins);    
    [Npro,binpro,idxpro]        = histcounts(kr_pro_all{i},kbins);    
    [Nlove,binlove,idxlove]     = histcounts(kr_love_all{i},kbins);
    
    % Considering beam response amplitude in wavenumber histogram
    if isempty(idxretro) == 0
        Aretro = accumarray(idxretro',a_retro_all{i});
        Aretro(end:numel(Nretro)) = 0;
        % Pick maximum using findpeaks
        Aretro_norm = Aretro./mean(Aretro); % normalise using mean
        [peak,idx,wpeak] = findpeaks(Aretro_norm,'NPeaks',1,'SortStr','descend','MinPeakHeight',SNR,'WidthReference','halfheight');
        if isempty(idx)==0
            kh_retro_max(i) = kbins(idx); % wavenumber of peak
            wmax_retro_k(i) = wpeak * kstep/2; % halfwidth of peak
            prmask(i) = i; % index of peak
        end
        k_retro_amp(i,:) = Aretro;
        k_retro_norm(i,:) = Aretro_norm;
        clear Aretro*
    end
    
    if isempty(idxpro) == 0
        Apro = accumarray(idxpro',a_pro_all{i});
        Apro(end:numel(Npro)) = 0;
        % Pick maximum using findpeaks
        Apro_norm = Apro./mean(Apro); % normalise using mean
        [peak,idx,wpeak] = findpeaks(Apro_norm,'NPeaks',1,'SortStr','descend','MinPeakHeight',SNR,'WidthReference','halfheight');
        if isempty(idx)==0
            kh_pro_max(i) = kbins(idx); % wavenumber of peak
            wmax_pro_k(i) = wpeak * kstep/2; % halfwidth of peak
            ppmask(i) = i; % index of peak
        end
        k_pro_amp(i,:) = Apro;
        k_pro_norm(i,:) = Apro_norm;
        clear Apro*
    end
    
    if isempty(idxlove) == 0
        Alove = accumarray(idxlove',a_love_all{i});
        Alove(end:numel(Nlove)) = 0;
        % Pick maximum using findpeaks
        Alove_norm = Alove./mean(Alove); % normalise using mean
        [peak,idx,wpeak] = findpeaks(Alove_norm,'NPeaks',1,'SortStr','descend','MinPeakHeight',SNR,'WidthReference','halfheight');
        if isempty(idx)==0
            kh_love_max(i) = kbins(idx); % wavenumber of peak
            wmax_love_k(i) = wpeak * kstep/2; % halfwidth of peak
            plmask(i) = i; % index of peak
        end
        k_love_amp(i,:) = Alove;
        k_love_norm(i,:) = Alove_norm;
        clear Alove
    end
    
end

%% Plot wavenumber histograms
    
% Retrograde Rayleigh

if plotkhisto
    
cmap = crameri('batlow');

if exist('prmask','var')==0
    warning('Retrograde Rayleigh waves: no peak detected above SNR!')
else

    figure('Color','white','Position',[2000 10 1000 500])
    % imagesc(freqs,kbins,kh_retro');
    imagesc(freqs,kbins,k_retro_amp');
    set(gca,'ydir','normal')
    xlabel('frequency in Hz')
    ylabel('wavenumber in 1/m')
    set(gca,'xtick',0.2:0.1:1.2,'Fontsize',18)
    colormap(cmap)
    cr = colorbar;
    % cr.Label.String = 'counts';
    cr.Label.String = 'cummulative amplitude';
    title('f-k histogram: retrograde Rayleigh waves')
    hold on
    errorbar(freqs(prmask > 0),kh_retro_max(prmask > 0),wmax_retro_k(prmask > 0),'w')
    % plot(freqs,kr_retro_mean,'w','Linewidth',2)
    % plot(freqs,kretro_max,'w--','Linewidth',2)
    
    legend('maximum')
    lhisto=legend('boxoff');
    lhisto.TextColor = [1 1 1];
    
    if savekhisto
        fname = sprintf('%s/khisto_retro_%s_%s_%s.png',...
            dir_out,maxflag,nw,indate);
        export_fig(fname,'-dpng');
end

end

% Prograde Rayleigh

if exist('ppmask','var')==0
    warning('Prograde Rayleigh waves: no peak detected above SNR!')
else
    
    figure('Color','white','Position',[2000 10 1000 500])
    imagesc(freqs,kbins,k_pro_amp');
    set(gca,'ydir','normal')
    xlabel('frequency in Hz')
    ylabel('wavenumber in 1/m')
    set(gca,'xtick',0.2:0.1:1.2,'Fontsize',18)
    colormap(cmap)
    cr = colorbar;
    % cr.Label.String = 'counts';
    cr.Label.String = 'cummulative amplitude';
    title('f-k histogram: prograde Rayleigh waves')
    hold on
    errorbar(freqs(ppmask > 0),kh_pro_max(ppmask > 0),wmax_pro_k(ppmask > 0),'w')
    % plot(freqs,kr_pro_mean,'w','Linewidth',2)
    % plot(freqs,kpro_max,'w--','Linewidth',2)
    legend('maximum')
    lhisto=legend('boxoff');
    lhisto.TextColor = [1 1 1];
    
    if savekhisto
        fname = sprintf('%s/khisto_pro_%s_%s_%s.png',...
            dir_out,maxflag,nw,indate);
        export_fig(fname,'-dpng');
    end

end

% Love

if exist('plmask','var')==0
    warning('Love waves: no peak detected above SNR!')
else

figure('Color','white','Position',[2000 10 1000 500])
imagesc(freqs,kbins,k_love_amp'); 
set(gca,'ydir','normal')
xlabel('frequency in Hz')
ylabel('wavenumber in 1/m')
set(gca,'xtick',0.2:0.1:1.2,'Fontsize',18)
colormap(cmap)
cr = colorbar;
% cr.Label.String = 'counts';
cr.Label.String = 'cummulative amplitude';
title('f-k histogram: Love waves')
hold on
errorbar(freqs(plmask > 0),kh_love_max(plmask > 0),wmax_love_k(plmask > 0),'w')
legend('maximum')
lhisto=legend('boxoff');
lhisto.TextColor = [1 1 1];

if savekhisto
    fname = sprintf('%s/khisto_love_%s_%s_%s.png',...
        dir_out,maxflag,nw,indate);
    export_fig(fname,'-dpng');
end

end 

end % if plotkhisto

%% Plot velocity dispersion curves

if plotvdisp
    
% Compute velocities from maxima in wavenumber histogram

% Retrograde Rayleigh
if exist('prmask','var')==1
vretro = freqs(prmask > 0)./kh_retro_max(prmask > 0);
wmax_retro_v = wmax_retro_k(prmask > 0)./kh_retro_max(prmask > 0) .* vretro;
end

% Prograde Rayleigh
if exist('ppmask','var')==1
vpro = freqs(ppmask > 0)./kh_pro_max(ppmask > 0);
wmax_pro_v = wmax_pro_k(ppmask > 0)./kh_pro_max(ppmask > 0) .* vpro;
end

% Love
if exist('plmask','var')==1
vlove = freqs(plmask > 0)./kh_love_max(plmask > 0);
wmax_love_v = wmax_love_k(plmask > 0)./kh_love_max(plmask > 0) .* vlove;
end

figure('Color','white','Position',[2000 10 1000 500])
legend('-DynamicLegend')
hold all
if exist('prmask','var')==1; e1 = errorbar(freqs(prmask>0),vretro,wmax_retro_v,...
        'LineWidth',2,'Color',batlowS(1,:),'DisplayName','Rayleigh (retro)'); end
if exist('ppmask','var')==1; e2 = errorbar(freqs(ppmask>0),vpro,wmax_pro_v,...
        'LineWidth',2,'Color',batlowS(2,:),'DisplayName','Rayleigh (pro)'); end
if exist('plmask','var')==1; e3 = errorbar(freqs(plmask>0),vlove,wmax_love_v,...
        'LineWidth',2,'Color',batlowS(3,:),'DisplayName','Love/SH'); end
title('Dispersion curves')
xlabel('frequency in Hz')
ylabel('velocity in m/s')
set(gca,'Fontsize',18)
box on
grid on
ylim([0 5000])
xlim([freqs(1) freqs(end)])

if savevdisp
    fname = sprintf('%s/disp_%s_%s_%s.png',...
        dir_out,maxflag,nw,indate);
    export_fig(fname,'-dpng');
end

end

%% Plot direction of arrival (azimuth)

if plotDOA

cmap = crameri('batlow');    

% Histograms
clear N*
abins = -175:5:180;
for i = 1:length(freqs)
    
    [Nretro,binretro,idxretro] = histcounts(kth_retro_all{i}*180/pi,abins);    
	[Npro,binpro,idxpro] = histcounts(kth_pro_all{i}*180/pi,abins);
	[Nlove,binlove,idxlove] = histcounts(kth_love_all{i}*180/pi,abins);
    [Nsv,binsv,idxsv] = histcounts(kth_sv_all{i}*180/pi,abins);
    [Np,binp,idxp] = histcounts(kth_p_all{i}*180/pi,abins);

    % Considering beam response amplitude in DOA histogram
    if isempty(idxretro) == 0
        Aretro = accumarray(idxretro',a_retro_all{i});
        Aretro(end:numel(Nretro)) = 0;
        azi_retro_amp(i,:) = Aretro;
        clear Aretro*
    end
    
    if isempty(idxpro) == 0
        Apro = accumarray(idxpro',a_pro_all{i});
        Apro(end:numel(Npro)) = 0;
        azi_pro_amp(i,:) = Apro;
        clear Apro*
    end
    
    if isempty(idxlove) == 0
        Alove = accumarray(idxlove',a_love_all{i});
        Alove(end:numel(Nlove)) = 0;
        azi_love_amp(i,:) = Alove;
        clear Alove*
    end
    
    if isempty(idxsv) == 0
        Asv = accumarray(idxsv',a_sv_all{i});
        Asv(end:numel(Nsv)) = 0;
        azi_sv_amp(i,:) = Asv;
        clear Asv*
    end
    
    if isempty(idxp) == 0
        Ap = accumarray(idxp',a_p_all{i});
        Ap(end:numel(Np)) = 0;
        azi_p_amp(i,:) = Ap;
        clear Ap*
    end
end

xmax = max(freqs)+0.2*max(freqs);

figure('Color','white','Position',[50 50 1400 700]);

% Retrograde Rayleigh
tldoa = tiledlayout(2,3);
nexttile
Nretroshift = circshift(azi_retro_amp,round(size(azi_retro_amp,2)/2),2);
%     Nretroshift = circshift(Nretro,round(size(Nretro,2)/2),2);
polarplot3d(Nretroshift,'RadialRange',[min(freqs) max(freqs)],'PlotType','surfn','PolarDirection','ccw','TickSpacing',15);
view(0,90)
colormap(cmap)
axis equal
box on
grid off
set(gca,'Fontsize',18)
title('Retrograde Rayleigh wave')
%     set(gca,'xtick',[],'ytick',[])
xlim([-xmax xmax])
ylim([-xmax xmax])
xlabel('frequency in Hz')
ylabel('frequency in Hz')
cb = colorbar;
cb.Label.String = 'cummulative amplitude';

% Prograde Rayleigh
nexttile
Nproshift = circshift(azi_pro_amp,round(size(Npro,2)/2),2);
polarplot3d(Nproshift,'RadialRange',[min(freqs) max(freqs)],'PlotType','surfn','PolarDirection','ccw','TickSpacing',15);
view(0,90)
colormap(cmap)
axis equal
box on
grid off
set(gca,'Fontsize',18)
title('Prograde Rayleigh wave')
%     set(gca,'xtick',[],'ytick',[])
xlim([-xmax xmax])
ylim([-xmax xmax])
xlabel('frequency in Hz')
ylabel('frequency in Hz')
cb = colorbar;
cb.Label.String = 'cummulative amplitude';

% SH/Love wave
nexttile
Nloveshift = circshift(azi_love_amp,round(size(Nlove,2)/2),2);
polarplot3d(Nloveshift,'RadialRange',[min(freqs) max(freqs)],'PlotType','surfn','PolarDirection','ccw','TickSpacing',15);
view(0,90)
colormap(cmap)
axis equal
box on
grid off
set(gca,'Fontsize',18)
title('SH/Love wave')
% set(gca,'xtick',[],'ytick',[])
xlim([-xmax xmax])
ylim([-xmax xmax])
xlabel('frequency in Hz')
ylabel('frequency in Hz')
cb = colorbar;
cb.Label.String = 'cummulative amplitude';

% P wave
nexttile
Npshift = circshift(azi_p_amp,round(size(Np,2)/2),2);
polarplot3d(Npshift,'RadialRange',[min(freqs) max(freqs)],'PlotType','surfn','PolarDirection','ccw','TickSpacing',15);
view(0,90)
colormap(cmap)
axis equal
box on
grid off
set(gca,'Fontsize',18)
title('P wave')
% set(gca,'xtick',[],'ytick',[])
xlim([-xmax xmax])
ylim([-xmax xmax])
xlabel('frequency in Hz')
ylabel('frequency in Hz')
cb = colorbar;
cb.Label.String = 'cummulative amplitude';

% SV wave
nexttile
Nsvshift = circshift(azi_sv_amp,round(size(Nsv,2)/2),2);
polarplot3d(Nsvshift,'RadialRange',[min(freqs) max(freqs)],'PlotType','surfn','PolarDirection','ccw','TickSpacing',15);
view(0,90)
colormap(cmap)
axis equal
box on
grid off
set(gca,'Fontsize',18)
title('SV wave')
% set(gca,'xtick',[],'ytick',[])
xlim([-xmax xmax])
ylim([-xmax xmax])
xlabel('frequency in Hz')
ylabel('frequency in Hz')
cb = colorbar;
cb.Label.String = 'cummulative amplitude';

if saveDOA
    fname = sprintf('%s/DOA_%s_%s_%s.png',...
        dir_out,maxflag,nw,indate);
    export_fig(fname,'-dpng');
end

end

%% EOF

