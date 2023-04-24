function [kr_max, kth_max, pola_max, pola_ind, a_max, wave_ind] = f_extrema22(P, Q, kr, kth, polstates, crit1, crit2, min_beam)

nwin = size(P,3);

kr_max = zeros(nwin,1);
kth_max = zeros(nwin,1);
pola_max = zeros(nwin,1,4);
a_max = zeros(nwin,1);
pola_ind = zeros(nwin,1);
wave_ind = zeros(nwin,1);

% Loop over time windows
for itime = 1:nwin
        
    P1 = P(:,:,itime);
    Q1 = Q(:,:,itime);
    
    pmean = mean(mean(P1));
    pstd = std(std(P1));
        
    % Find strongest peaks in the spectrum
    [ax,idx] = extrema2(P1); % get amplitudes and indices of global extrema
    
    % make sure peak is larger than mean + 3 std (background "noise")
    if ax(1) < pmean+3*pstd
        continue
    end
    
    % First criterion: minimum amplitude
    if strcmp(crit1, 'MIN')
        % Discard maximum when amplitude < min_beam (ratio of first
        % maximum)
        % --> see Riahi et al. 2013 p. 5343
        ax_max = ax(1);
        idx = idx(ax >= min_beam * ax_max & ax > pmean+3*pstd);
        ax = ax(ax >= 0.5 * ax_max & ax > pmean+3*pstd);
    elseif strcmp(crit1, 'NOMIN')
        % Use maxima regardless of amplitude;
        % only make sure it's above noise level
        idx = idx(ax > pmean+3*pstd);
        ax = ax(ax > pmean+3*pstd);
    end
    
    [i,j] = ind2sub(size(P1),idx);
    
    % Correct for double count at edges (0 and 360 degree)
    i1 = find(i==1);
    i2 = find(i==length(kth),1);
    if ~isempty(i1) && ~isempty(i2)        
        if P1(i(i1),j(i1)) > P1(i(i2),j(i2))
            j = j(i~=length(kth));
            i = i(i~=length(kth));
        else
            j = j(i~=1);
            i = i(i~=1);
        end
        P1_help = P1(:);
        idx = sub2ind(size(P1),i,j);
        ax = P1_help(idx);%/P1_help(idx_help(1));
    end
    
    if isempty(i)
        disp(itime)
    end
    
    % Second criterion: maximum number of events
    if strcmp(crit2, 'NOMAX')
        % Use all events detected
        i = i;
        j = j;
    elseif strcmp(crit2,'MAX1')
        % Use first maximum only
        i = i(1);
        j = j(1);
    elseif strcmp(crit2, 'MAX3')
        % Use first three maxima
        i = i(1:min(3,length(idx)));
        j = j(1:min(3,length(idx)));
    elseif strcmp(crit2, 'MAX5')
        % Use first five maxima
        i = i(1:min(5,length(idx)));
        j = j(1:min(5,length(idx)));
    end
    ni = length(i);
    
    kr_max(itime,1:ni) = kr(j);%*1000; % GIVES OUTPUT IN 1/KM!
    kth_max(itime,1:ni) = kth(i);
    a_max(itime,1:ni) = ax(1:ni);
    
    for c = 1:ni
        
        pola_help = polstates(Q1(i(c),j(c)),:);
        
        pola_max(itime,c,:) = pola_help;
        pola_ind(itime,c) = Q1(i(c),j(c));
        
        if pola_help(3) == 0 && pola_help(4) == pi
            % P-wave
            wave_ind(itime,c) = 0;
        elseif pola_help(3) == 2 && pola_help(4) == pi
            % SV-wave
            wave_ind(itime,c) = 1;
        elseif pola_help(3) == 2 && pola_help(4) == pi/2
            % Love wave
            wave_ind(itime,c) = 2;
        elseif pola_help(3) > 0 && pola_help(3) < 2 && pola_help(4) == 0
            % Retrograde Rayleigh wave
            wave_ind(itime,c) = 3;
        elseif pola_help(3) > 0 && pola_help(3) < 2 && pola_help(4) == pi
            % Prograde Rayleigh wave
            wave_ind(itime,c) = 4;
        end
        
    end
      
end % time windows

end % function