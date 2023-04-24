function [mindist, maxdist, imin, jmin, imax, jmax] = f_compminmaxdist(coords)
% Compute minimum and maximum interstation distance within an array
% required: coords (n x 2) double

mindist = 10^5;
maxdist = 0;
for i = 1:size(coords,1)
    for j = 1:size(coords,1)
        if i ~= j
            dist = norm(coords(i,:)-coords(j,:));
            if dist < mindist
                mindist = dist;
                imin = i;
                jmin = j;
            end
            if dist > maxdist
                maxdist = dist;
                imax = i;
                jmax = j;
            end
        end
        
    end
end
