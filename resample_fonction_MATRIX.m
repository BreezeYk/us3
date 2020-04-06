function [signal_surech,t_surech] = resample_fonction_MATRIX(signal,Fs,F)
% Cette fonction permet de rééchantillonner toutes les colonnes d'une 
% matrice pour faire passer les signaux de Fs à F sur la même durée 
% temporelle. Une interpolation de type linéaire est utilisée. 
%
% Entrees : - signal = signal à rééchantillonner, chaque colonne = 1 signal
%           - Fs = frequence d'échantillonnage du signal 
%           - F = fréquence d'échantillonnage voulue
% Sorties : - signal_surech = signal rééchantillonné
%           - t_surech = temps associe à signal_surech
%
% Eloise Chassaing - 12/7/2018

% signal
[duree,nb_rec] = size(signal); 

% temps associé au signal en entrée :
t_fs=[0:1:duree-1]/Fs;

% temps souhaité :
t_surech=[0:Fs/F:duree-1]/Fs; 

% on interpole :
for i=1:nb_rec
    signal_surech(:,i) = interp1(t_fs,signal(:,i),t_surech,'linear');
end

end

