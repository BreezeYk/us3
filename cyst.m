% clc;clear all;close all;

% definition parametres
% load('../data/paramModeleDirecte/paramModeleDirecte.mat');
function [pt_positions amp] = cyst(N)

%      N=15e3; % definition nombre des diffuseurs
%     N = 5.5195e3 ;
    z_min=0.01;   % definition de la position z minimale
    z_max=0.02;   % definition de la position z maximale

    x_min=-0.005;       % definition de la position x minimale
    x_max=0.005;        % definition de la position x maximale

    xc_1=0;
    zc_1=0.015;

    xc_2=0;
    zc_2=0.030;

    xc_3=0.0045;
    zc_3=0.035;

    rayon_123=0.002;

    profondeur=0.04;
    eppaisseur=0.005;

    pt_positions=[rand(N,1) zeros(N,1) rand(N,1)];

    pt_positions=bsxfun(@plus,[x_min 0 z_min],bsxfun(@times,[(x_max-x_min) 0 (z_max-z_min)],pt_positions));


    coeffs_1=((pt_positions(:,1)-xc_1).^2+(pt_positions(:,3)-zc_1).^2)<rayon_123^2;
    % coeffs_2=((pt_positions(:,1)-xc_2).^2+(pt_positions(:,3)-zc_2).^2)<rayon_123^2;
    % coeffs_3=((pt_positions(:,1)-xc_3).^2+(pt_positions(:,3)-zc_3).^2)<rayon_123^2;

    coeffs_4=((pt_positions(:,3)-profondeur)<=eppaisseur)&((pt_positions(:,3)-profondeur)>0);

    amp=(rand(N,1)-0.5)*2;
    amp(coeffs_1)=0;
    % amp(coeffs_2)=2*sign(amp(coeffs_2))+amp(coeffs_2);
    % amp(coeffs_3)=1*sign(amp(coeffs_3))+amp(coeffs_3);

    amp=amp./max(abs(amp));

    amp(coeffs_4)=min(amp)+(max(amp)-min(amp))*(pt_positions(coeffs_4,1)-x_min)/(x_max-x_min);
    display_scatterers(pt_positions,amp,[min(pt_positions(:,1)) max(pt_positions(:,1))], [min(pt_positions(:,3)) max(pt_positions(:,3))]);

end

    

    % save('../data/phantoms/large_full_medium_fancy.mat','pt_positions','amp');