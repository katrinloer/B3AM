function mat = smoothmat(mat,smkernel)
% 
% Smooth a matrix along its rows and/or columns. Essentially a conv2
% operation, but with appropriate cropping to keep the size of the data
% matrix the same. Smoothing window (1xn, nx1, nxm) must be supplied.
% 
%   mat = nri_smoothmat(mat,smkernel)
% 
% Typically:
% smkernel = ones(1,nsmo)/nsmo;
% OR
% hwin = hanning(nsmo);
% smkernel = hwin * hwin';
% smkernel = smkernel/sum(smkernel(:));
% OR similar...
% 


[nsmoC,nsmoR] = size(smkernel);
if min(nsmoC,nsmoR)~=0
    if nsmoC==1 % row smooth
        mat = conv2(mat,smkernel);  % Order of arguments of conv2 does not matter
        mat = mat(:,1+(nsmoR-1)/2:end-(nsmoR-1)/2);
    elseif nsmoR==1
        mat = conv2(mat,smkernel);  
        mat = mat(1+(nsmoC-1)/2:end-(nsmoC-1)/2,:);
    else
        mat = conv2(mat,smkernel);  
        mat = mat(1+(nsmoC-1)/2:end-(nsmoC-1)/2 , 1+(nsmoR-1)/2:end-(nsmoR-1)/2);
    end
end



end