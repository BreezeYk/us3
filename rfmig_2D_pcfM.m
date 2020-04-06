function [migSIG1,varargout] = rfmig_2D_pcf(SIG,x,z,param,varargin)

%   rfmig_2D_pcf   Migration of Raw signals (with Phase Coherence Factor in option)
%   Edited version of rfmig_2D() (from US_toolbox from GIT)
%   migSIG1 could be:
%   - migrated signals (param.do_summation=0;)
%   - RF image (param.do_summation=1;)
%   SIG could be:
%   - RF raw signals
%   - IQ signals (computed with rf2iq() from US_toolbox from GIT)
%   param is a structure that contains all the parameter values required
%   for the migration (see below for more details).

%   * Dynamic aperture can be used (param.dyn_rx_apert=1);
%   In this case, you can get edges of dynamic aperture (in number of
%   elemnt) for each pixel in ind_bornes (size Nz, Nx, 2): 
%   [migSIG1,~,ind_bornes] = rfmig_2D_pcf(raw,x,z,param)
%   Required when using MVbeamforming2 for spatial smoothing.
%
%   * Adaptive Phase Coherence Factor can be computed (param.do_pcf =1)
%   Then, multiply the DAS image obtain with the weighting "image" wPCF
%   [migSIG1,wPCF] = rfmig_2D_pcf(raw,x,z,param)
%
%   * Delayed Phases to later compute PCF can be returned(param.do_pcf =1)
%   [migSIG1,~,~,phase] = rfmig_2D_pcf(raw,x,z,param)
%
%%% Faire un meshgrid des vecteurs x et z 

%   The signals - typically RF signals - in SIG must be acquired using a
%   full wave (planar, circular, parabolic...) configuration as used in
%   ultrafast ultrasound imaging. Each column corresponds to a single RF
%   signal over time (as acquired by a single transducer element).
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%   1) PARAM.fs: sample frequency (in Hz, REQUIRED)
%   2) PARAM.pitch: pitch of the transducer (in m, REQUIRED)
%   3) PARAM.xn: virtual source position of x direction(in m, REQUIRED)
%      PARAM.yn: virtual source position of y direction(in m, REQUIRED)
%      PARAM.zn: virtual source position of z direction(in m, REQUIRED)
%   4) PARAM.c0: longitudinal velocity (in m/s, default = 1540 m/s)
%   5) PARAM.t0: start time for reception (in s, default = 0 s)
%   6) PARAM.dyn_rx_apert : '1' if you want dynamic aperture size (default), '0' 
%   7) PARAM.do_summation : '1' if you want Delay and sum (default), '0' if you just
%   want delayed signals
%   8) PARAM.do_pcf : '1' if you want to get the PCF weighting mask 
%	9) param.xm / ym / zm = position de la sonde
%
% Last modification : M. Polichetti (10/12/2019)

assert(ndims(SIG)<=3,['SIG must be a matrix whose each column ',...
    'corresponds to an RF signal acquired by a single element']);
varargout{1}=0;
varargout{2}=0;
%-- Check input parameters
sizx = size(x);
if nargin==4
    assert(isequal(size(z),sizx),'X and Z must be of same size.')
end

if (~isfield(param,'c0'))
    param.c0 = 1540;
end

if ~isfield(param,'fs')
    error('A sampling frequency (PARAM.fs) is required.')
end

if ~isfield(param,'t0')
    param.t0 = 0; % acquisition start time in s
end

if ~isfield(param,'Pitch_x') % in m
    if isfield(param,'dx')
        param.Pitch_x = param.dx; % param.dx was used in the old version
    elseif isfield(param,'Pitch')
        param.Pitch_x = param.Pitch; % param.dx was used in the old version
    elseif isfield(param,'pitch')
        param.Pitch_x = param.pitch; % param.dx was used in the old version
    else
        error('A pitch value (PARAM.Pitch_x) is required.')
    end
end

if ~isfield(param,'xn')
    param.xn = 0;
end

if ~isfield(param,'yn')
    param.yn = 0;
end

if ~isfield(param,'zn')
    param.zn = -10;
end

if ~isfield(param,'F')
    param.F = 3;
end

if ~isfield(param,'RXangle')
    param.RXangle = 0;
end

x1 = param.xm;  % Element positions
y1 = param.ym;
z1 = param.zm; 

c = param.c0;   % propagation velocity

if ~isfield(param, 'type')
    param.type = 'TX';
end

if ~isfield(param, 'dyn_rx_apert')
    param.dyn_rx_apert = 0;
end

if ~isfield(param, 'do_summation')
    param.do_summation = 1;
end

if ~isfield(param, 'do_pcf')
    param.do_pcf = 0;
end

if ~isfield(param, 'do_oversampling')
    param.do_oversampling = 0;
    warning('Please, check if param.fs>8*param.f0, to avoid reconstruction artefacts. Or, give me param.do_oversampling=1')
end

if (param.do_oversampling==1)&&(isreal(SIG))
    if(~isfield(param, 'f0'))
        warning('param.f0 is missing to do oversampling. Default value for oversampling factor is 2.')
        param.f0=param.fs/4;
    end
%     disp('Coucou')
    if (8*param.f0>param.fs)
        previous_fs = param.fs;     % save parameter
        n_ov = ceil(8*param.f0/param.fs); % Ovsersampling factor param.fs/previous_fs
        [SIG, param.fs] = oversample_fourier(SIG,previous_fs,n_ov);
    end
end

%-- Migration 
if (param.do_summation)
    migSIG = zeros(size(x));
else
    migSIG = zeros([numel(x1) numel(x)]);
end

%-- Migration 
if (param.dyn_rx_apert==1)
    ind_borne = zeros([numel(x),2]);
else 
    ind_borne = ones(numel(x),2);
    ind_borne(:,2) = length(x1)*ind_borne(:,2); 
    if nargout==3
        varargout{2} = reshape(ind_borne,[size(x),2]);
    end
end

% compute the minimal distance from the probe to the virtual source
dmin = min( sqrt( (x1(:)-param.xn).^2 + (y1(:)-param.yn).^2 + (z1(:) -param.zn).^2 ) );
if (param.zn<0)
    dmin = -dmin;
end

%% Delay And Sum

if param.do_pcf == 1
    mask_pcf = zeros([length(x1) numel(x) 2]);      % Element, pixels, (phi, phi_c)
    if (isreal(SIG))  % Select Phase on analytic signals
        phi = angle(hilbert(SIG));     % Phases des rawdata
        
    else % If IQ, re-introduce phase on IQ
        t = ((0:size(SIG,1)-1)).'/param.fs;  %-- Time vector
        phi = angle(SIG.*exp(2*1i*pi*param.f0*t));
    end
    % Init. variables
    PHI = zeros(numel(x),1);
    PHI_c = PHI;
    gam0 = 1;
    sig0 = pi/3^0.5;
end

for k = 1:numel(x1)     % Loop on the elements
    
    %-- Compute distances for all pixels
    dTX = sqrt((x(:)-param.xn).^2 + (z(:)-param.zn).^2) + dmin; % transmit distance
    dRX = sqrt((x1(k)-x(:)).^2 + (z(:).^2)); % receipt distance
    
    %-- Convert distances into delays
    if (strcmp(param.type, 'RX'))
        tau = (dRX) / c;
    else
        tau = (dTX + dRX) / c;
    end
    
    %-- Convert delays into samples
    idxt = (tau  - param.t0)*param.fs +1;
%     idxt = (tau - param.td - param.t0)*param.fs +1;
    I = idxt<1 | idxt>(size(SIG,1)-1);
    idxt(I) = 1; % arbitrary index, will be soon rejected
     
    idx  = idxt;       % Not rounded number of samples to interpolate later
    idxf = floor(idx); % rounded number of samples 
    IDX  = repmat(idx, [1 1 size(SIG,3)]);
    
    %-- Recover delayed samples with a linear interpolation 
    TEMP = SIG(idxf,k,:).*(idxf+1-IDX) + SIG(idxf+1,k,:).*(IDX-idxf);
    
    if (~isreal(TEMP)) % if IQ signals
        TEMP = TEMP.*exp(2*1i*pi*param.f0*tau);
    end
    
    TEMP(I,:) = 0;
    
    %-- Reshape data
    if(param.dyn_rx_apert==0)
        TEMP = reshape(TEMP, [size(x) size(SIG,3)]); % Just Reshape
    else
        %-- Reshape and perform dynamic aperture with a Constant F-Number 
        if (param.RXangle==0)
            ind = abs(x-x1(k)) < z/param.F/2;
        else
            ind = abs(x-x1(k)+tan(param.RXangle)*z) < z/cos(param.RXangle)/param.F/2;
        end
        TEMP = repmat(ind,[1 1 size(SIG,3)]) .* reshape(TEMP, [size(x) size(SIG,3)]); 
        
        % Store boundaries of dynamic aperture
        if (param.dyn_rx_apert==1)
            ind_borne((ind(:)==1)&(ind_borne(:,1)==0),1) = k;
            ind_borne(:,2) = ind_borne(:,2) + ind(:);
            if k == numel(x1)
                ind_borne(:,2)=ind_borne(:,2)+ind_borne(:,1)-1; 
                if nargout==3
                    varargout{2} = reshape(ind_borne,[size(x),2]);
                end
            end
        end
        
    end
    
    if param.do_pcf==1
        % Select and interpolate the Phases
        PHI = phi(idxf,k,:); %phi(idxf,k,:).*(idxf+1-IDX) + phi(idxf+1,k,:).*(IDX-idxf);
        PHI_c(PHI>0) = PHI(PHI>0)-pi; % Shift the phases
        PHI_c(PHI<0) = PHI(PHI<0)+pi; % Shift the phases
        mask_pcf(k,:,1) = PHI;   % Store the phases
        mask_pcf(k,:,2) = PHI_c; % Store the shifted phases
    end
    
    % Store or Sum the Focused line
    if (param.do_summation)
        migSIG = migSIG + TEMP ; % Sum
    else
        migSIG(k,:) = TEMP(:);   % Store
    end
    
end

    if (param.do_pcf==1)
        % Compute the weights for PCF for all pixels
        var_phi   = std(mask_pcf(:,:,1),1);     % Variance of the phase of the aperture for each pixel
        var_phi_c = std(mask_pcf(:,:,2),1);     % Variance of shifted phase of the aperture for each pixel
        [var_phi_f,nphi] = min([var_phi(:),var_phi_c(:)]'); % Keep the the lower for each pixel
        PCF_w = reshape(max(zeros(numel(var_phi_f'),1), 1-gam0/sig0.*var_phi_f'),size(x)); % Compute de the PCF for each pixel
        varargout{1}=PCF_w;
                       
         % To return phases and not just PCF
        if (nargout>=4)
            % Select phases (corrected or not) used to compute PCF
            NN = repmat(nphi,[length(x1),1]);
            mask_pcf_new = ((NN==1).*mask_pcf(:,:,1)) + ((NN==2).*mask_pcf(:,:,2));
            varargout{3} = reshape(mask_pcf_new.',[size(x),length(x1)]); % [z,x,elem]
        end
    end

    % Reshape data before return
    if (param.do_summation)
       migSIG1=reshape(migSIG, size(x));      
    else
        migSIG1=reshape(migSIG, [length(x1) size(x)]);
        migSIG1=permute(migSIG1,[2 3 1]);

    end

end

