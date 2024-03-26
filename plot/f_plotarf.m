function [] = f_plotarf(coords,kmin,kmax,kres,ares,scaleflag,save_figs,figpath)

n = size(coords,1);

%% Plot array
figure('Color','white','Position',[10 400 500 500]);
plot(coords(:,1),coords(:,2),'kv','MarkerFaceColor','k','MarkerSize',12)
box on
grid on
axis equal
xlabel('Easting in m')
ylabel('Northing in m')
set(gca,'Fontsize',20)
title('Array geometry')

if save_figs
    figname = [figpath, '/01_network.png'];
    print(figname,'-dpng')
end

%% Array response vector
kr = linspace(0,kmax,kres)';
kth = (-(180-ares):ares:180)*pi/180;
k = kron(kr, [cos(kth(:)) sin(kth(:))] );

ak = 1/sqrt(n) * exp(1i*2*pi*coords*k');

%% Vertical incidence wavefield
kr_s = 0; % wave coming from below
kth_s = kth;
k_s = kron(kr_s, [cos(kth_s(:)) sin(kth_s(:))]);
s = exp(1i*2*pi*coords*k_s'); % wavefield of different k_s
s = mean(s,2); % average over waves with different k_s (coming from different azimuths)

%% Array response:
% R(k) = a(k)* (s s*) a(k), where (s s*) = S (spectral density matrix)
% see Riahi et al. (2013) eq. 3
R = conj(ak)' * s;
R = R .* conj(R);

R = reshape(R,size(kth,2),size(kr,1));

if strcmp(scaleflag,'norm')
    R = R / max(max(R)); % normalize
elseif strcmp(scaleflag,'db')
    R = mag2db(R); % convert magnitude to decibel
end

%% Plot ARF

if kmax < 10^-3
    kmaxp = kmax * 1000;
    kminp = kmin * 1000;
    krp = kr * 1000;
    labelunit = 'km^{-1}';
else
    kmaxp = kmax;
    kminp = kmin;
    krp = kr;
    labelunit = 'm^{-1}';
end

% Response (polar)
figure('Color','white','Position',[110 400 500 500]);
Rcirc=circshift(R',round(size(R',2)/2),2);
polarplot3d(Rcirc,'PlotType','surfn','RadialRange',[0 kmaxp],'PolarDirection','ccw',...
    'RadLabels',0,'RadLabelLocation',{90, 'top'},'AxisLocation','off','Fontsize',12,'RadLabelColor','w');
% polarplot3d(Rcirc,'PlotType','surfn','RadialRange',[0 kmax*1000],'PolarDirection','ccw',...
%     'RadLabels',4,'RadLabelLocation',{90, 'top'},'AxisLocation','off','Fontsize',12);
view(0,90);
axis equal
xlim([-kmaxp kmaxp])
ylim([-kmaxp kmaxp])
set(gca,'Fontsize',20)%'XTick',-0.3:0.1:0.3)
box on
xlabel(sprintf('k_x in %s',labelunit))
ylabel(sprintf('k_y in %s',labelunit))
caxis([0 1])
title({'ARF'})
% colormap('jet')
cmap = crameri('batlow');
colormap(cmap);
cbar = colorbar;
if strcmp(scaleflag,'norm')
    cbar.Label.String = 'normalized beam power';
elseif strcmp(scaleflag,'db')
    cbar.Label.String = 'beam power in decibel';
end
cbar.FontSize = 20;
cbar.Label.FontSize = 20;

if save_figs
    figname = [figpath, '/02_ARF.png'];
    print(figname,'-dpng')
end

%% Plot ARF cross-sections (see Wathelet et al. 2008)

figure('Color','white','Position',[210 400 500 500]);
p1 = plot(krp,R,'Color',[0.5 0.5 0.5]);
hold on; 
l1 = line([0 kmaxp],[0.5 0.5],'Color','k','Linewidth',2,'LineStyle',':');
l2 = line([kminp kminp],[0 1],'Color','k','Linewidth',2);
xlabel(sprintf('wavenumber in %s',labelunit))
ylabel('normalized beampower')
xlim([0 kmaxp])
box on
grid on
set(gca,'Fontsize',20)
title('ARF cross-sections')
legend([l1, l2], 'half height', 'min. wavenumber')

if save_figs
    figname = [figpath, '/03_ARFx.png'];
    print(figname,'-dpng')
end

end