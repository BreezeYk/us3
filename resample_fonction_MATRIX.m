function [signal_surech,t_surech] = resample_fonction_MATRIX(signal,Fs,F)
% Cette fonction permet de r��chantillonner toutes les colonnes d'une 
% matrice pour faire passer les signaux de Fs � F sur la m�me dur�e 
% temporelle. Une interpolation de type lin�aire est utilis�e. 
%
% Entrees : - signal = signal � r��chantillonner, chaque colonne = 1 signal
%           - Fs = frequence d'�chantillonnage du signal 
%           - F = fr�quence d'�chantillonnage voulue
% Sorties : - signal_surech = signal r��chantillonn�
%           - t_surech = temps associe � signal_surech
%
% Eloise Chassaing - 12/7/2018

% signal
[duree,nb_rec] = size(signal); 

% temps associ� au signal en entr�e :
t_fs=[0:1:duree-1]/Fs;

% temps souhait� :
t_surech=[0:Fs/F:duree-1]/Fs; 

% on interpole :
for i=1:nb_rec
    signal_surech(:,i) = interp1(t_fs,signal(:,i),t_surech,'linear');
end

end

