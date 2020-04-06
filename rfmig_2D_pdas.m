function [im_pdas, varargout] = rfmig_2D_pdas(SIG,X,Z,h,p, varargin)
% Function to reconstruct 'im_pdas', RF-Image using pth root compression
% Signals are migrated using rfmig_2D_pcf()
    % SIG raw RF signals (depth,elements) to be migrated, or already
    % migrated (then h.do_mig = 0)
    % x lateral coordinates for migration (meshgrid)
    % z depth coordinates for migration (meshgrid)
    % p range for pth root compression
    % h settings structure for rfmig_2D_pcf
    % wn : (optional) normalized cutoff frequencies. Default: [0.3 1.7]*h.f0/(h.fs/2)
% /!\ A previous oversampling of SIG (fs>8f0) may be required
% Maxime Polichetti, Jan.,2018

h.do_pcf = 0;
h.do_summation = 0;

if ~isfield(h,'do_mig')
    h.do_mig = 1;
end

%-- Migrate raw signals
    if h.do_mig == 1
        migSIG = rfmig_2D_pcf(SIG, X, Z, h);
    else
        migSIG = SIG;
    end
    
%-- pth rooth compression
    migSIG_temp = sign(migSIG).*(abs(migSIG).^(1/p));
  
%-- Sum on the element dimension
    im_pdas_temp = sum(migSIG_temp,3);
    
%-- p-Powering to recover original dimensonnality
    im_pdas_temp = sign(im_pdas_temp).*(abs(im_pdas_temp).^p);
     
% Band pass Filter to reject induced artificial harmonics
if (nargin>5)
    wn = varargin{1};
else
    wn = [0.3 1.7]*h.f0/(h.fs/2);
end

    if wn ==0
        im_pdas = im_pdas_temp;
        
    else

        [b,a] = butter(11,wn(2));
        im_pdas = filtfilt(b,a,im_pdas_temp);

        if (wn(1)~=0)
            [b,a] = butter(11,wn(1),'high');
            im_pdas = filtfilt(b,a,im_pdas);
        end
    end
    
if nargout==2 % Do simple DAS
    varargout{1} = sum(migSIG,3); % im_das
end

end

