function [IQ,Fc] = rf2iq(RF,Fs,Fc,Wn)

%RF2IQ   I/Q demodulation of RF data
%   IQ = RF2IQ(RF,Fs,Fc) demodulates the radiofrequency (RF) bandpass
%   signals and returns the Inphase/Quadrature (I/Q) components. IQ is a
%   complex whose real (imaginary) part contains the inphase (quadrature)
%   component.
%       1) Fs is the sampling frequency of the RF signals (in Hz),
%       2) Fc represents the center frequency (in Hz).
%
%   IQ = RF2IQ(RF,Fs) or IQ = RF2IQ(RF,Fs,[],...) determines the carrier
%   frequency automatically.
%     -- !WARNING! Fc must be given if the RF signal is undersampled --
%
%   RF2IQ uses a downmixing process followed by low-pass filtering. The
%   low-pass filter is determined by its normalized cut-off frequency (Wn):
%   IQ = RF2IQ(RF,Fs,Fc,Wn) uses Wn as the normalized cut-off frequency.
%   By default, Wn = 2*Fc/Fs. The lower bound for Wn is B/Fs, where B is
%   the bandwidth of the bandpass RF signal.
%
%   [IQ,Fc] = RF2IQ(...) also returns the carrier frequency (in Hz).
%
%   Notes on Undersampling (sub-Nyquist sampling)
%   ----------------------
%   If the RF signal has been undersampled, the carrier frequency Fc
%   (before undersampling) must be given. A warning message appears if
%   harmful aliasing is suspected.
%
%   Notes:
%   -----
%   RF2IQ treats the data along the first non-singleton dimension as
%   vectors, i.e. RF2IQ demodulates along columns for 2-D and 3-D RF data.
%   Use IQ2RF to recover the RF signals.
%
%   Method:
%   ------
%   RF2IQ multiplies RF by a phasor of frequency Fc (down-mixing) and
%   applies a fifth-order Butterworth lowpass filter using FILTFILT:
%       IQ = RF.*exp(-1i*2*pi*Fc*t);
%       [b,a] = butter(5,2*Fc/Fs);
%       IQ = filtfilt(b,a,IQ)*2;
%
%   References:
%   ----------
%   1) J Kirkhorn, Introduction to IQ-demodulation of RF-data, 1999. 
%   <a
%   href="matlab:web('http://folk.ntnu.no/htorp/Undervisning/TTK10/IQdemodulation.pdf')">PDF download</a>
%   2) SA Tretter, Bandpass signal representation, 1999. 
%   <a
%   href="matlab:web('http://www.ece.umd.edu/class/enee429w.F99/bandpass.pdf')">PDF download</a>
%   3) M Mishali and YC Eldar, Sub-Nyquist sampling, 2011. 
%   <a
%   href="matlab:web('http://webee.technion.ac.il/Sites/People/YoninaEldar/journals/06021873-123.pdf')">PDF download</a>
%
%   Example #1
%   ----------
%   % Load an RF signal sampled at 20 MHz
%   load RFsignal@20MHz.mat
%   % I/Q demodulation
%   IQsignal = rf2iq(RFsignal,20e6);
%   % RF signal and its envelope
%   plot(RFsignal), hold on
%   plot(abs(IQsignal),'Linewidth',1.5), hold off
%   legend({'RF signal','I/Q amplitude'})
%
%   Example #2: Demodulation of an undersampled RF signal
%   ----------
%   % Load an RF signal sampled at 20 MHz
%   % (Center frequency = 5 MHz / Bandwidth = 2 MHz)
%   load RFsignal@20MHz.mat
%   % I/Q demodulation of the original RF signal
%   Fs = 20e6; Fc = 5e6;
%   IQ = rf2iq(RFsignal,Fs,Fc);
%   % Create an undersampled RF signal (sampling at Fs/5 = 4 MHz)
%   rf = RFsignal(1:5:end);
%   % I/Q demodulation of the undersampled RF signal
%   Fs = 4e6; B = 2e6; Wn = B/Fs;
%   iq = rf2iq(rf,Fs,Fc,Wn);
%   % Display the results
%   plot(abs(IQ(1:5:end)),'Linewidth',1.5), hold on
%   plot(abs(iq),'r'), hold off
%   title('I/Q amplitude (5 MHz array)')
%   legend({'sampled @ 20 MHz','undersampled @ 4 MHz'})
%
%   See also IQ2RF.
%
%   -- Damien Garcia -- 2012/01, last update: 2016/07
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>


%-- Check input arguments
narginchk(2,4);
assert(isreal(RF),'RF must contain real signals.')
assert(isscalar(Fs),'Fs (sampling frequency in Hz) must be a scalar.')
if nargin>2
    assert(isempty(Fc) || isscalar(Fc),'Fc must be [] or a scalar.')
end

%-- Convert to column vector (if RF is a row vector)
if isrow(RF)
    RF = RF(:);
    wasrow = true;
else
    wasrow = false;
end

%-- Time vector
nl = size(RF,1);
t = (0:nl-1)'/Fs;

%-- Seek the carrier frequency (if required)
if nargin<3 || isempty(Fc)
    % Keep a maximum of 100 scanlines
    Nc = size(RF,2);
    if Nc<100, idx = 1:Nc; else idx = randperm(Nc); idx = idx(1:100); end
    % Power Spectrum
    P = sum(abs(fft(RF(:,idx))).^2,2);
    P = P(1:floor(nl/2)+1);
    % Carrier frequency
    idx = sum((0:floor(nl/2))'.*P)/sum(P);
    Fc = (idx-1)*Fs/nl;
end

%-- Normalized cut-off frequency
if nargin<4 || isempty(Wn)
    Wn = 2*Fc/Fs-eps;
end
assert(isscalar(Wn),'Wn (normalized cutoff frequency) must be a scalar.')
assert(Wn>0 && Wn<1,'Wn (normalized cutoff frequency) must be within the interval of (0,1).')

%-- Down-mixing of the RF signals
IQ = bsxfun(@times,double(RF),exp(-1i*2*pi*Fc*t));

%-- Low-pass filter
[b,a] = butter(5,Wn);
IQ = filtfilt(b,a,IQ)*2; % factor 2: to preserve the envelope amplitude

%-- Recover the initial size (if was a vector row)
if wasrow, IQ = IQ.'; end

%-- Display a warning message if harmless aliasing is present
B = Fs*Wn; % bandwidth
fL = Fc-B/2; fH = Fc+B/2; % lower and higher frequencies of the bandpass signal
n = floor(fH/(fH-fL));
harmlessAliasing = any(2*fH./(1:n)<=Fs & Fs<=2*fL./(0:n-1));
if ~harmlessAliasing
    warning('MUT:harmfulAliasing',...
        'Harmful aliasing is present: the aliases are not mutually exclusive!')
end
