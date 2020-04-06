% Function to oversample a signal by zero-padding its spectrum.
% Extension of 'interpft.m' Matlab function, to non integer oversampling factor
%
% To use it :
%   [s_out, fs_out]        = oversample_fourier( s_in, fs_in, n_ov);
%   [s_out, fs_out, t_out] = oversample_fourier( s_in, fs_in, n_ov, t_in);
% with:
%   - s_in is the initial signal  
%     (if s_in matrix, first dimension is the temporal one)
%   - fs_in is the initial sampling frequency (Hz)
%     (fs_in = 1, if digital frequency)
%   - n_ov is the approximative oversampling factor fs_out/fs_in
%   - s_out is the oversampled version s_in
%   - fs_out is the sampling frequency of s_out 
%
% Note that the first sample for both s_in and s_out is computed for the
% same initial time. The new temporal axis t_out can be automatically
% computed adding the initial temporal axis t_in as input.
% 
% Example :
%         % Create Signal
%             fs = 18;
%             t  = 0 : 1/fs : 1 - 1/fs;
%             s  = cos(2*pi*3*t);
%         % Upsample it
%             n_ov = 4.5;
%             [s_new, fs_new, t_new] = oversample_fourier( s, fs, n_ov, t);
%         % Display Signals
%             figure, 
%             stem(t, s, 'b'), hold on
%             stem(t_new, s_new, 'r--')
%             xlabel('Time (s)')
%             title('Signals')
%   
% M. Polichetti : maxime.polichetti@creatis.univ-lyon1.fr

function [ s_out, fs_out, varargout ] = oversample_fourier( s_in, fs_in, n_ov, varargin)

% Works with column vector
    if(isrow(s_in)==1)
        s_in=s_in';
    end

% Number of samples
    N_in    = size(s_in,1);   % initial number of samples
    N_out   = ceil(n_ov*N_in);   % number of samples after oversampling

% Constraint on constant frequency spacing for both spectrum
    df      = fs_in/N_in;   % frequency spacing
    fs_out  = df * (N_out); % exact oversampling frequency       

% Compute the FFT of initial signal
    TFs_in   = fft(s_in);    

% Construct what would be the FFT of the oversampled signal
    % Positive frequencies until fs/2 : TFs values
    % Positive frequencies above fs/2 : 0 (zero-padding)
    % Negative frequencies : arbitrary set to 0 (cf ifft(,symmetric))
    TFs_out = zeros(N_out,size(s_in,2),size(s_in,3));
    TFs_out(1:floor(N_in/2)+1,:,:) = TFs_in(1:floor(N_in/2)+1,:,:);
    
% Go back to the time domain
    s_out = ifft(n_ov*TFs_out,'symmetric');
    
% Avoid discontinuities (update, July 2018)
    [b,a] = butter(5,(fs_in/2)/(fs_out/2));
    s_out = filtfilt(b,a,double(s_out));
%     warning('Check filter - oversample_fourier')

% Return new temporal vector if asked
    if (nargin==4) && (nargout==3)
        t_in = varargin{1};    % time vector for s
        t_out = t_in(1) + (0:1/fs_out:(N_out-1)/fs_out);
        varargout{1} = t_out;
    end
end