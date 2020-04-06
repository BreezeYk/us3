function [ s_out ] = add_noise( s_in, SNR_dB, paramField , sig_type )
% Function to add noise to data s_in [Time, Element], with a signal to noise
% ratio SNR_dB in dB, which is automatically adapted considering the
% data_type: "Stable" for Narrowband signal on 1 bin, else, s_in is
% considered as broadband signal.
% varargin{1} is a seed to fix random number generator
% varargin{2} is the frequency sampling for s_in required for plot

    
% Make the difference between narrowband and broadband signals
% and compute the number of frequency bins on which to add noise
if strcmp('Stable',sig_type)
    N_bin_signal = 1;
    N_bin_bruit  = round(size(s_in,1)/2);
    N_snr = N_bin_signal/N_bin_bruit;
else
    N_bin_signal = round(size(s_in,1)/2);
    N_bin_bruit  = round(size(s_in,1)/2);
    N_snr = N_bin_signal/N_bin_bruit;
end
    
% Compute the mean power of data s_in
    Psignal     = mean(mean(s_in.^2,1));
    
% Convert SNR_dB into a linear ratio, to deduce the corresponding std for
% the white gaussian noise
    SNR_lin     = 10^(SNR_dB/10);
    std_noise   = sqrt(Psignal/SNR_lin/N_snr);
    fn = 0.5*paramField.fs;
    fd = 0.1e6;
    fu = 20e6;
    [B, A] = fir1(50,[fd fu]/fn);
    noise       = std_noise*randn(size(s_in));
    noise_limited = filter(B, A, noise);
    

% Add noise
    s_out    = s_in + noise_limited;    
    
% Display
%    if(nargin >=5)
 %       fs  = varargin{2};
        
  %      s   = s_in;
   %     b   = noise;
    %    TFs = abs(1/length(s)*fft(s));
    %    TFb = abs(1/length(b)*fft(b));
     %   axf = 0 : fs/length(TFs) : (fs-fs/length(TFs));
      %  
       % figure
       % plot(axf, 20*log10(mean(TFb,2)), 'r','Linewidth',1.5)
       % hold on
       % plot(axf, 20*log10(mean(TFs,2)), 'b','Linewidth',1.5)
       % xlim([0 fs/2])
       % legend('Bruit','Signal')
       % xlabel('Frequence')
  %  end
    
end

