% EXAMPLE
% clear all, close all, clc

% addpath(genpath('../'))
if exist( 'dataset' ) ~= 1
  %  load Beamforming_Dataset.mat
end
% Raw data (1PW with 0° angle)
%     SIG = dataset.data(:,:,3);
%SIG = [];

%SIG = h_hatt( 1:length(h_hatt)/paramEm.K , : ) ; 

%load('h_est_2scat_inf_dB_K4.mat') ;
a = size( squeeze( h_target(1,:,:) ) , 1 ) - size( h_hatt( 1:length(h_hatt)/paramEm.K , : ) , 1 ) ;
 SIG1 = cat( 1 , squeeze(h_est(1,:,:)) , zeros(a,64)) ;
% SIG1 = cat( 1 , squeeze(mf(:,:,1))' , zeros(a,64)) ;

%load('h_est_2scat_5_dB_K4.mat') ;
a = size( squeeze( h_target(1,:,:) ) , 1 ) - size( h_hatt( 1:length(h_hatt)/paramEm.K , : ) , 1 ) ;
 SIGb = cat( 1 , squeeze(h_est(1,:,:)) , zeros(a,64)) ;
 
 SIGr = squeeze( h_target(1,:,:) )  ;
%SIG = ( h_est' )  ;
% for i = 1: size( signauxRF1.rfSignals , 2)
  %  for i = 1: 16
    
   % SIG = [SIG  signauxRF1.rfSignals{i} ] ;
    
%        SIG = [SIG a(:,i)] ;
    
%end


% Plug the header
    h.fs    = paramField.fs;
    h.pitch = abs(paramField.probe_specs(8,2) - paramField.probe_specs(8,1));
%     h.pitch = paramField.pitch ;
    h.c0    = paramField.c;
    h.f0    = paramField.f0;
    h.N_cycle = 2; % cf Challenge PICMUS spec.
%     h.t0    = signauxRF2.offsetTime{1} - h.N_cycle/h.f0/2;
   % h.t0    = signauxRF1.offsetTime{1} - h.N_cycle/h.f0/2;
     h.t0    = os - h.N_cycle/h.f0/2;
    % Source (here 1PW)
    h.xn = 0;    h.yn = 0;    h.zn = -10;
    % Elements (here linear array)
%     h.xm = dataset.probe_geometry(:,1);    h.ym = 0*h.xm;    h.zm = 0*h.xm;
    h.xm = paramField.probe_specs(8,:);    h.ym = 0*h.xm;    h.zm = 0*h.xm;
    % What you want to do
    h.dyn_rx_apert = 1;    h.F = 3;
    h.do_summation = 1;
    h.do_pcf= 0;
    h.do_oversampling = 1;
%    h.td = temps_decalage_final(1) ;
    
% Reconstruction axis
x = linspace(h.xm(1),h.xm(end),ceil(length(h.xm)/2)*2+1);
z = (0:size(SIG)-1)*h.c0/h.fs/2;
idx = find( z < 0.06 ) ;


[X,Z] = meshgrid(x,z);

% Reconstruct RF image with DAS and p-DAS
p=0; % range for pth root compression
%[im_pdas , im_das1 , W , A , B ] = rfmig_2D_pdasM( SIG ,X,Z,h,p , 0 );

[im_das_1] = rfmig_2D(SIG1,X,Z,h,p);
[im_das_B] = rfmig_2D(SIGb,X,Z,h,p);
[im_das_r] = rfmig_2D(SIGr,X,Z,h,p);


im_das_pcf = im_das1.*W;
[b,a]      = butter(11,[0.5 1.5]*h.f0/(h.fs/2));
im_das_pcf = filtfilt(b,a,double( im_das_pcf));

%% Display


imr = sum( im_das_r , 1 ) ;
imb = sum( im_das_B , 1 ) ;
im1 = sum( im_das_1 , 1 ) ;
dyn = 60; 
figure
% subplot(121),
subplot(131);    imagesc(x,z,rf2bmode(imr(1,:,:))) ; title('Ref') ;
subplot(132);    imagesc(x,z,rf2bmode(imb(1,:,:))) ; title('Noisy') ;
subplot(133);    imagesc(x,z,rf2bmode(im1(1,:,:))) ; title('Non noisy') ;

% imagesc(x,z(idx),rf2bmode(im_pdas(idx,:))) ;
% imagesc(x,z(idx),rf2bmode(im_das_pcf(idx,:))) ;
% imagesc(x,z(idx),rf2bmode( im_das(8,idx,:))) ;

title('DAS'), axis image
% subplot(122), imagesc(x,z,rf2log(im_pdas,dyn)), title('p-DAS'), axis image
colormap gray
%%
im1 = sum( im_das_1 , 1 ) ;
imB = sum( im_das_2 , 1 ) ;
imNoise = imB - im1 ;

figure;
 imagesc(x,z,rf2bmode(imNoise(1,:,:))) ; axis image ; colormap gray ;
 
s = rf2bmode(im(1,:,:)) ;