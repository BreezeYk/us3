function [chirp_compresse,chirp_compresse2] = compression_chirp_filtre_wiener3(rf_signal, excitation,alpha)
% ATTENTION ! CETTE FONCTION NECESSITE LA CONNAISSANCE DU PARAMETRE ALPHA !
%
% Fonction : Comprimer les signaux rf (sous forme de matrice) par un 
% filtre de Wiener. 
% Méthode utilisée : multiplication fréquentielle. 
% 
% Entrées : - rf_signal : matrice de signaux RF de taille
%           L0*Nelements_reception, correspond aux chirps réceptionnés lors
%           d'une émission
%           - excitation : signal d'excitation (chirp long temporellement),
%           vecteur ligne (taille 1*n)DOUBLEMENT CONVOLUE AVEC LA REPONSE
%           - alpha : paramètre utilisateur
%           IMPULSIONNELLE DES TRANSDUCTEURS avec la fonction conv
% Sorties : - chirp_compresse : chirp comprimé temporellement par filtre 
%           de Wiener (taille : L0*Nelements_reception) de manière à 
%			ce que le pic de compression coïncide avec le début du chirp 
%			- chirp_compresse2 : chirp comprimé par la méthode
%			fréquentielle, de manière à ce que le signal comprimé 
%			commence avant l'arrivée du chirp (informations plus complètes 
%			mais temps des données à modifier en conséquence, voir le script)
%
% 
% Eloise Chassaing - 23/06/2018

if size(excitation,2) > 1
    excitation = excitation'; 
end

[L0, nb_recepteurs] = size(rf_signal); 
L = max(L0,length(excitation)); 
if L < 2^11
    L = 2^11;  
end
L = 2^nextpow2(L);
S = fft(rf_signal, L); % spectre du signal reçu bruité

%% signal d'excitation convolué 2 fois avec la réponse impulsionnelle du
% transducteur
% on construit Ah = transfo de Fourier de l'excitation doublement convoluée
% if length(rep_impuls) == 1
%     Ah = (fft(rep_impuls,L).*fft(rep_impuls,L))'.*fft(excitation',L);
% elseif ~isempty(rep_impuls)
%     Ah = fft(rep_impuls',L).*fft(rep_impuls',L).*fft(excitation',L);
% else
%     Ah = fft(excitation',L);
% end
Ah=fft(excitation,L); 
% ah = ifft(Ah);
Ah = repmat(Ah,1,nb_recepteurs); 

% cas où il n'y a pas la double convolution : la compression marche (pas de
% décalage temporel)
[ah_comp] = compression_chirp_filtre_adapte(excitation, excitation');
% figure,plot(ah_comp)

% on cherche le décalage temporel dû à la double convolution, qui sera
% soustrait à nos résultats de filtrage
% [ah_comp] = compression_chirp_filtre_adapte(ah, excitation);
% env_ah_comp = abs(hilbert(ah_comp)); 
% 
% if ~isempty(rep_impuls)
%     [~,max_pic_ah_comp] = max(env_ah_comp);
% else
%     max_pic_ah_comp = 1; 
% end

%%
% alpha = 50; %% faire une boucle pour trouver le meilleur : attention si 
% trop petit, filtre va être instable (car se rapproche du filtre inverse) 
% donc va amplifier les petits trucs (bruit)

%% on construit le filtre de Wiener : 
H = conj(Ah)./ ((abs(Ah)).^2 + alpha); 
% freq_h = (0:length(H)-1)*fs/length(H)/10^6; 
% figure, plot(freq_h,abs(H)),title('Réponse fréquentielle du filtre (avec convolution)'); 

%% on filtre nos données :
filtrage_freq = S.*H; 
filtrage_temporel = real(ifft(filtrage_freq(:,:))); 
% figure,plot(filtrage_temporel(:,1))

% on les décale temporellement pour prendre en compte la double convolution
% [c,~]=find(excitation ~= 0); 
% c=min(c); 
% chirp_compresse = [(zeros(c-1,nb_recepteurs)) ; filtrage_temporel(1:L0,:) ];% (zeros(max_pic_ah_comp-1,nb_recepteurs)) ; 

chirp_compresse = [(zeros(0,nb_recepteurs)) ; filtrage_temporel(1:L0,:) ];% (zeros(max_pic_ah_comp-1,nb_recepteurs)) ; 

% 
% figure, 
% plot(chirp_compresse(:,1))
%/max(chirp_compresse(:,1)),'r'), hold on, plot(rf_signal(:,1)/max(rf_signal(:,1)),'b'),title('Résultat de la compression avec décalage');


%% si on ne veut pas couper les chirps : 
% on garde tout le signal mais on renvoie là où il faut couper
% 
l = size(filtrage_temporel,1); 
chirp_compresse2 = [ filtrage_temporel(l-floor(length(excitation))+1:end,:) ; filtrage_temporel(1:l-floor(length(excitation)),:)];
% [ filtrage_temporel(floor(L/2)+1:end,:) ; filtrage_temporel(1:floor(L/2),:)];


% % donc la troncature va être : 
% temps = (length(excitation))+1; %floor(L/2); 
% chirp_compresse2 = chirp_compresse2(temps:end,:) ;
% figure, plot(chirp_compresse2(:,1))
end