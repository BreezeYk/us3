% This function convert a radio-frequency image into a log-compress image.
% To use it :
%      Im_log = rf2log(rf_in)
%      Im_log = rf2log(rf_in, dB)
% where 
%   - rf_in is the input RF image
%   - Im_log is the output log-compress image
%   - dB is an optional argument to impose the dnamic of the log-compressed
%   image. If not given, 40dB is used

function Im_log = rf2log(rf_in, varargin)
if (nargin==2)
    Im_log = bmode2log(rf2bmode(rf_in), varargin{1});
else
    Im_log = bmode2log(rf2bmode(rf_in));
end