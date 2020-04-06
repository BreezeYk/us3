% This function convert a radio-freuqncy image into a bmode image.
% To use it :
%      Bmode = rf2bmode(rf_in)
% where 
%   - rf_in is the input RF image
%   - Bmode is the output bmode image
function rf_out = rf2bmode(rf_in)
rf_in = double( squeeze(rf_in) );
if (length(size(rf_in))<3)
    rf_out = abs(hilbert(rf_in));
elseif (length(size(rf_in))==3)
    rf_out = 0*rf_in;
    for i=1:size(rf_in,3)
        rf_out(:,:,i) = abs(hilbert(rf_in(:,:,i)));
    end
else
    error('The dimension of the input image is incorrect');
end