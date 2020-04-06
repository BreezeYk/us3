function [chirp_compresse_adapte3,chirp_compresse_adapte] = compression_chirp_filtre_adapte(rf_signal, excitation)
% Fonction : Comprimer les signaux rf (sous forme de matrice) par un 
% filtre adapté. 
% 3 méthodes possibles sont utilisées : convolution temporelle,
% multiplication fréquentielle et intercorrélation. 
% 
% Entrées : - rf_signal : matrice de signaux RF de taille
%           L*Nelements_reception, correspond aux chirps réceptionnés lors
%           d'une émission
%           - excitation : signal d'excitation (chirp long temporellement),
%           vecteur ligne (taille 1*n) DOUBLEMENT CONVOLUE
% Sorties : - chirp_compresse_adapte3 : chirp comprimé par la méthode
%           fréquentielle (taille L*Nelements_reception), de manière à 
%			ce que le pic de compression coïncide avec le début du chirp 
%			- chirp_compresse_adapte : chirp comprimé par la méthode
%			fréquentielle, de manière à ce que le signal comprimé 
%			commence avant l'arrivée du chirp (informations plus complètes 
%			mais temps des données à modifier en conséquence, voir le script)
% Autres sorties possibles : 
%           - chirp_compresse_adapte2 : chirp comprimé temporellement par
%           intercorrelation (taille L*Nelements_reception)
%           - chirp_compresse_adapte1 : chirp comprimé temporellement par
%           convolution (taille L*Nelements_reception)
% 
% Eloise Chassaing - 11/04/2018

% on met le signal d'excitation 
if size(excitation,1)>1
    excitation = excitation'; 
end

[L,nb_recepteurs] = size(rf_signal); 


% filtre adapté = signal d'excitation retourné temporellement
filtre_adapte = flip(excitation); % retournement temporel


%{
%%% Méthodes 1 et 2 (convolution temporelle et intercorrélation)
for i = 1: nb_recepteurs
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1ere methode : convolution : 
    tmp = conv(rf_signal(:,i),filtre_adapte);
%     figure, plot(tmp); 
    chirp_compresse_adapte1(:,i) = tmp(length(excitation):end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2ème méthode : intercorrélation
    [tmp2, ~] = xcorr(rf_signal(:,i),excitation); % la fonction rajoute des zéros si les 2 n'ont pas la meme longueur
    % figure, plot(tmp2) ;
    L0 = max(L,length(excitation)); 
    chirp_compresse_adapte2(:,i) = tmp2(length(excitation):end);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3ème méthode : en fréquentiel : 
% % comme les signaux n'ont pas forcément la meme longueur, on force à avoir
% % une longueur suffisante
L1 = max(size(rf_signal,1),length(excitation)); %2^12; %
if L1<1024
    L2=2*1024; 
else
    L2=2*L1;
end
X = fft(rf_signal,L2);  
Y = fft(filtre_adapte',L2) ;
Y2 = repmat(Y, 1,nb_recepteurs);
tmp3 = ((real(ifft(X.*Y2))));
chirp_compresse_adapte3 = tmp3(length(excitation):end,:); 

if size(chirp_compresse_adapte3,1)>=L
    chirp_compresse_adapte3=chirp_compresse_adapte3(1:L,:);
else
    chirp_compresse_adapte3=[chirp_compresse_adapte3 ; zeros(L-size(chirp_compresse_adapte3,1),nb_recepteurs)];
end

%%%% figure avec les 3 méthodes
% figure, plot(chirp_compresse_adapte3(:,i),'*r'); hold on, plot(chirp_compresse_adapte1(:,i),'b');hold on, plot(chirp_compresse_adapte2(:,i),'v');

    
%% si on ne veut pas couper les chirps : 
% on garde tout le signal mais on renvoie là où il faut couper
% renvoyer les tmp et dire où il faudra considérer le temps de début
% égal à length(excitation)
if L+length(excitation)-1<= size(tmp3,1)
    chirp_compresse_adapte=tmp3(1:L+length(excitation)-1,:);
else
    chirp_compresse_adapte=tmp3(1:L,:);
end
% donc la troncature va être : 
% chirp_compresse_adapte=chirp_compresse_adapte(length(excitation):end,:);
    
end