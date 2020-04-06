function [rfBeamformed, x_grille, z_grille, stockage_temps, num0,numm0,amplitude_point_suivi] = reconstruire_image_STA_codes_recepteurs(recordedData,paramImage, paramEmission, paramField, temps_decalage)
% Fonction :- Reconstruire les images basse r�solution � partir des signaux 
% RF r�ceptionn�s lors d'une ouverture synth�tique mono-�l�ment (un �l�ment
% �met et tous re�oivent). L'utilisateur peut pr�ciser les r�cepteurs qu'il
% veut utiliser lors de la reconstruction. 
%           - Visualiser les temps utilis�s, pour un point du fantome,
% dans le calcul de l'interpolation. 
%           - Si l'utilisateur pr�cise un point qui doit �tre au centre de
%           l'image, la fonction cr�e une grille centr�e sur ce point, avec
%           les r�solutions d�termin�es (param�tre nomm�
%           paramImage.diffuseur_reso)
% N.B : - Si le milieu imag� est moins profond que la grille demand�e, les
% param�tres de formation de l'image sont modifi�s pour ne traiter que les 
% profondeurs r�ellement enregistr�es ! Cela �vite la bande noire
% observ�e s'il n'y a pas de donn�es � interpoler. 
%       - Type d'interpolation : lin�aire
%       - Code issu de celui de Denis Bujoreanu
% 
%
% Entr�es : - recordedData : structure contenant les signaux r�ceptionn�s
%           avec 2 champs : rfSignals et offsetTime. 
%           Taille =  L*Nelements_reception
%           - paramImage : structure contenant les param�tres des images
%           basse r�solution � cr�er (x_min, x_max, resolution_laterale, 
%           prof_min_milieu, prof_max_milieu, resolution_axiale) et le
%           point du fantome a suivre lors de l'interpolation 
%           (coordonnees en m�tre dans position_a_suivre). Elle peut 
%           contenir aussi le champ diffuseur_reso (coordonn�es d'un point
%           qu'on veut au centre des images basse r�solution) (FACULTATIF)
%           - paramEmission : structure contenant les param�tres des
%           �missions des signaux (N_STA, sources_emettrices,
%           position_elements, recepteurs_actifs)
% recepteurs_actifs
%           - paramField : structure contenant les param�tres de Field II
%           (c et fs)
%           - temps_decalage : temps du maximum de l'enveloppe du signal
%           d'excitation doublement convolu�e (SI N'A PAS ETE PRIS EN COMPTE
%           DANS offsetTime de recordedData) (en seconde)
%
% Sorties : - rfBeamformed : cube (de taille x*z*N_STA)
%           contenant chacune des images basse r�solution reconstruites 
%           - x_grille : vecteur ligne de taille 1*x, contenant les
%           abscisses des points utilis�s dans la grille des images basse
%           r�solution
%           - z_grille : vecteur colonne de taille z*1, contenant les
%           profondeurs des points utilis�s dans la grille des images basse
%           r�solution
%           - stockage_temps : matrice de taille 
%           Nelements_reception*Nemissions qui contient sur chaque
%           colonne les temps utilis�s pour interpoler l'amplitude du point 
%           de la grille le plus proche de (x0,z0)(pour chaque image BR)
%           (1 colonne = 1 �mission = 1 image BR // 1 ligne = temps 
%           utilis�s pour 1 r�cepteur)
%           - num0 : indice de position dans x_grille du point (x0,z0) suivi
%           - numm0 : indice de position dans z_grille du point (x0,z0) suivi
%           - amplitude_point_suivi : suivi de l'amplitude du point (x0,z0)
%           au cours des reconstructions des images basse r�solution
%
%   Eloise Chassaing - 11/04/2018

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du minimum / maximum de la profondeur des signaux recus +
% ajout de z�ros � nos donn�es recues si la profondeur r�ceptionn�e
% est plus petite que celle definie dans la grille (via prof_max_milieu)

% 1) r�duire la profondeur de la grille si on n'a pas les donn�es
% correspondantes



stock_profondeur = []; 
stock_debut = []; 
for source_index = paramEmission.sources_emettrices % pour chaque �mission
    % profondeur maximale dans les donn�es
    data_depth = size(recordedData.rfSignals{1,source_index},1)*paramField.c/(2*paramField.fs)+...
                                    recordedData.offsetTime{1,source_index}*paramField.c/2; 
    % profondeur minimale dans les donn�es
    min_depth =  recordedData.offsetTime{1,source_index}*paramField.c/(2);
    
    if data_depth < paramImage.prof_max_milieu % si nos donn�es ne vont pas aussi loin que la profondeur de la grille
        stock_profondeur = [stock_profondeur data_depth]; % on stocke la valeur de la profondeur maximale des donn�es
    end
    
    if min_depth > paramImage.prof_min_milieu % si nos donn�es commencent plus tard que ce qu'on voudrait pour la grille
        stock_debut = [stock_debut min_depth]; % on stocke cette profondeur minimale enregistr�e
    end
end

% si on a des donn�es qui ne vont pas aussi loin que la grille demand�e, on
% r�duit la profondeur de la grille
if ~isempty(stock_profondeur)
    min_profondeur = min(stock_profondeur); 
    paramImage.prof_max_milieu = min_profondeur; 
end 

% 2) calculer la taille minimale des signaux RF
%{
N_data2 = Inf; % longueur la plus courte des signaux recus apr�s ajout de 0 
for source_index = 1:paramEmission.N_STA % :paramField.N_elements_reception
    data_depth = size(recordedData2.rfSignals{source_index},1)*paramField.c/(2*paramField.fs)+...
                                    recordedData2.offsetTime{source_index}*paramField.c/2;             
    if data_depth < paramImage.prof_max_milieu % on va rajouter des 0 apr�s nos donn�es recues
        recordedData2.rfSignals{source_index}=[recordedData2.rfSignals{source_index};...
                        zeros(round(abs(data_depth-paramImage.prof_max_milieu)*2*paramField.fs/paramField.c), paramField.N_elements_reception)];
    end
    N_data2 = min([N_data2 size(recordedData2.rfSignals{source_index},1)]);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D�finition de la grille de points utilis�s (points qu'on va rechercher
% dans nos donn�es recues)
if isempty(paramImage.diffuseur_reso) % si on ne veut pas une grille centr�e sur un point
    x_grille = paramImage.x_min:paramImage.resolution_laterale:paramImage.x_max ; 
    z_grille = paramImage.prof_min_milieu:paramImage.resolution_axiale:paramImage.prof_max_milieu ;
else % si on veut une grille centr�e sur un point
    % v�rification 
    if paramImage.diffuseur_reso(1,1) < paramImage.x_min || paramImage.diffuseur_reso(1,1) > paramImage.x_max ...
            || paramImage.diffuseur_reso(1,3) < paramImage.prof_min_milieu || paramImage.diffuseur_reso(1,3) > paramImage.prof_max_milieu 
        error('revoir les coordonn�es du diffuseur dans paramImage.diffuseur_reso'); 
    end
    
    % calcul des grilles centr�es sur le diffuseur qu'on veut au centre de
    % l'image
    x_grille_droite = paramImage.diffuseur_reso(1,1):paramImage.resolution_laterale:paramImage.diffuseur_reso(1,1)+abs((paramImage.x_max-paramImage.x_min))/2; 
    x_grille_gauche = paramImage.diffuseur_reso(1,1):-paramImage.resolution_laterale:paramImage.diffuseur_reso(1,1)-abs((paramImage.x_max-paramImage.x_min))/2; 
    x_grille=[fliplr(x_grille_gauche(2:end)) x_grille_droite];
    
    z_grille_droite = paramImage.diffuseur_reso(1,3):paramImage.resolution_axiale:paramImage.diffuseur_reso(1,3)+abs((paramImage.prof_max_milieu-paramImage.prof_min_milieu))/2; 
    z_grille_gauche = paramImage.diffuseur_reso(1,3):-paramImage.resolution_axiale:paramImage.diffuseur_reso(1,3)-abs((paramImage.prof_max_milieu-paramImage.prof_min_milieu))/2; 
    z_grille=[fliplr(z_grille_gauche(2:end)) z_grille_droite];
end

M = length(x_grille); % nombre de pixels en x
N = length(z_grille); % nombre de pixels en z

imageGridX = repmat(x_grille,N,1);    % on calcule les positions en x
imageGridZ = repmat(z_grille',1,M);   % on calcule les positions en z

% toutes ces �tapes remplacent le meshgrid cr�� par la fonction : 
% [X, Z]=meshgrid(x_discret,z_discret); % mais celle-ci a un d�calage
% d'indices donc pas bon � utiliser ! 

x = imageGridX(:); % concat�ner les colonnes du meshgrid dans un vecteur
z = imageGridZ(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrouver le point d'interet x0, z0 dans la grille 
% on retrouve le num�ro de cellule (ligne/colonne) du point d'int�r�t
% d�fini dans paramImage.position_a_suivre
x0 = paramImage.position_a_suivre(1,1);
z0 = paramImage.position_a_suivre(1,3); 

num0 = find(imageGridX(1,:) ==x0); 
numm0 = find(imageGridZ(:,1) ==z0); 
if isempty(num0)
    tmp = imageGridX(1,:)-repmat(x0,1,size(imageGridX,2));
    tmp(tmp <0) = 999999; 
    [~,num0] = min(tmp);
    if num0>2 &&norm(imageGridX(1,num0)-x0,2) >= norm(imageGridX(1,num0-1)-x0,2)
        num0 = num0 -1 ; 
    end
end
if isempty(numm0)
    tmp = imageGridZ(:,1)-repmat(z0,size(imageGridZ,1),1);
    tmp(tmp <0) = 999999; 
    [~,numm0] = min(tmp);
    if numm0>=2 && norm(imageGridZ(numm0,1)-z0,2) >= norm(imageGridZ(numm0-1,1)-z0,2)
        numm0 = numm0 -1 ; 
    end
end

% trouvons la position de (x0,z0) dans x et z 
index0 = (num0-1) * size(imageGridZ,1) + numm0; 

% initialisation de la matrice qui va contenir les temps utilis�s pour
% interpoler le point (x0,z0)
stockage_temps = zeros( length(paramEmission.recepteurs_actifs) , paramEmission.N_STA); 
amplitude_point_suivi = []; % vecteur qui va contenir les amplitudes du 
% point suivi dans les images basse r�solution reconstruites  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de l'interpolation des temps et des amplitudes

f = waitbar(0,'Avancement du calcul');

tir = 1; % compteur du nb de tirs

for source_index = paramEmission.sources_emettrices  % pour chaque element ayant tir�
    
    waitbar(tir/paramEmission.N_STA,f, sprintf('Avancement de %d sur %d',[tir paramEmission.N_STA]));
    
    %%%%%%%%% TEMPS D'ALLER-RETOUR DES POINTS DE LA GRILLE %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Matrice de distance aller-retour
    mat_distance_aller = bsxfun(@plus,bsxfun(@minus,paramEmission.position_elements(source_index),x).^2,z.^2);  
    % attention, ici on peut avoir � prendre en compte des temps de retard
    % d'�mission (si source virtuelle ou onde plane)
    mat_distance_retour = bsxfun(@plus,bsxfun(@minus,paramEmission.position_elements,x).^2,z.^2);
    mat_distance_parcourue = bsxfun(@plus,sqrt(mat_distance_aller),sqrt(mat_distance_retour)); 
    
    % Conversion matrice de distance en matrice de temps
    mat_temps_parcours = mat_distance_parcourue / paramField.c ; 
    % ce sont nos temps de r�f�rence aller-retour pour les pts de la grille
    % pour l'�metteur consid�r�    
    
    % Prise en compte dans ces temps de r�f�rence du temps d'offset 
    % pr�sent dans nos donn�es recordedData2 et du temps de mont�e de l'enveloppe.            
    mat_temps_parcours = bsxfun(@minus, mat_temps_parcours, recordedData.offsetTime{1,source_index}-temps_decalage );
    % en effet, il faut faire au temps t : t-toffset+tdecalage. 
    
    % Conversion en nb d'�chantillons (plus pratique pour l'interpolation
    % lin�aire dans la suite)
    temps_parcours_ech = mat_temps_parcours*paramField.fs + 1; 
    
    % correspond aux temps d'aller-retour pour tous les points de la grille
    % pour chaque r�cepteur de la sonde (contient les temps vrais, absolus)
    % mais attention, ce ne sont pas forc�ment des nombres entiers
    % d'�chantillons !! alors que les temps des donn�es correspondent
    % forc�ment � des nbs entiers car les signaux RF bruts sont acquis
    % toutes les Ts secondes. 
    
    %%%%%%%%% REORGANISATION DES TEMPS DES POINTS DE LA GRILLE %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prise en compte de l'effet dent de scie pour la future interpolation
    % : besoin d'avoir un temps d'aller-retour qui augmente au fil de la
    % colonne
    N_data2 = Inf; 
    N_data2 = min([N_data2 size(recordedData.rfSignals{1,source_index},1)]);
    temps_parcours_ech = temps_parcours_ech(:,paramEmission.recepteurs_actifs); % on ne selectionne que nos recepteurs d'interet
    mat_temps_parcours_dent_scie_comp = bsxfun(@plus,temps_parcours_ech,(0:length(paramEmission.recepteurs_actifs)-1)*(N_data2));
    % permet d'avoir les n� de lignes correspondant quand on va tout mettre
    % sous forme d'un vecteur colonne dans l'interpolation entre les
    % donn�es RF enregistr�es et les temps des donn�es par rapport aux
    % temps de la grille

    %%%%%%%%%%%%%%%%%% TEMPS DANS LES DONNEES %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % L�, on va traiter les temps obtenus pour voir ceux qui sont 
    % effectivement r�alistes par rapport � nos donn�es RF re�ues
    % Ainis, on ne conserve que les temps >=1Ts et inf�rieur au temps des 
    % donn�es compl�t�es par des 0
    time_idx_in_data_in_bounds2 = temps_parcours_ech>= 1 & temps_parcours_ech<=(N_data2-1);
    % C'est une matrice qui ne contient que des 0 et des 1 pour dire quels
    % temps sont dans les donn�es
    
    % D�termination du vecteur temps des points de la grille (not� t ou vect_t)
    vect_t = mat_temps_parcours_dent_scie_comp(time_idx_in_data_in_bounds2);
    % c'est un vecteur 1D
        
    % On cherche dans notre matrice de temps de parcours en �chantillon 
    % avec prise en compte des dents de scie, les temps des donn�es, not� 
    % t_0 ou vect_t0. 
    vect_t0 = floor(vect_t); % floor -> correspond � de vrais �chantillons 
    % acquis avec les signaux RF (en fait des multiples de Ts). 
    
%     vect_t0
    % Explication de l'interpolation lineaire : 
    % Avec l'�galit� des coeff directeurs : 
    % (y(t)-y(t0))/(t-t0) = (y(t1)-y(t0))/(t1-t0)
    % Or on choisit t1 et t0 des temps successifs, donc en nombre
    % d'�chantillon, cela signifie que t1 - t0 = 1. 
    % On va donc faire le calcul : 
    % y(t) = y0(t0-t+1) - y1(t0-t) o� y(t) repr�sente l'amplitude du signal
    % RF recherch�e ! 
    
    % Calcul de t_0 - t 
    t0moinst = vect_t0 - vect_t;
                                       
    % Calcul de y(t) = y0(t0-t+1) - y1(t0-t)
%     data = recordedData.rfSignals{1,source_index}(:,paramEmission.recepteurs_actifs); 
     data = recordedData.rfSignals{1,source_index}(:,1); 
%      size(data)
%      size( vect_t0 )
%      data = recordedData.rfSignals{1,source_index}; 
    data_interp2 = data(vect_t0).*(t0moinst+1)-...
                            data(vect_t0+1).*t0moinst;
    % Stockage de l'interpolation lin�aire
    data_interp_els2 = zeros(numel(x),length(paramEmission.recepteurs_actifs));
    data_interp_els2(time_idx_in_data_in_bounds2) = data_interp2;
    
%     data_interp2 = recordedData.rfSignals{1,source_index}(vect_t0).*(t0moinst+1)-...
%                             recordedData.rfSignals{1,source_index}(vect_t0+1).*t0moinst;
%     % Stockage de l'interpolation lin�aire
%     data_interp_els2 = zeros(numel(x),paramField.N_elements_reception);
%     data_interp_els2(time_idx_in_data_in_bounds2) = data_interp2;
%     disp(['Valeur amplitude point suivi : ' num2str(data_interp_els2(numm0 +(num0-1)*length(z_grille),:))]);
%     disp(['Valeur amplitude point suivi somm�e : ' num2str(sum(data_interp_els2(numm0 +(num0-1)*length(z_grille),:))) ]); 
    
    % on remet en forme pour trouver notre grille de points 
    rfBeamformed(:,:,tir) = reshape(sum(data_interp_els2,2),N,M);
     
    amplitude_point_suivi = [amplitude_point_suivi rfBeamformed(numm0,num0,tir)];
    
    %%%%%% Visualisation de l'interpolation sur le diffuseur_suivi %%%%%%%
    % visualisation via croix
%     temps_point0 = zeros(numel(x),paramField.N_elements_reception);
    % cette matrice est de taille B*Nrecepteurs et va contenir les temps utilises 
    % pour trouver l'amplitude du point0 (par interpolation), selon chaque r�cepteur. 
%     vect_t2 = mat_temps_parcours(time_idx_in_data_in_bounds2); % on         
    % retrouve les temps de parcours sous forme de nb d'�chantillons,
    % en ne gardant que ceux des donn�es. puis on compl�te la matrice de
    % temps
%     temps_point0(time_idx_in_data_in_bounds2) = vect_t2; 
    % on convertit en temps (en seconde) : 
%     temps_point0 = temps_point0(index0,:)+recordedData.offsetTime{1,source_index}-temps_decalage;
    temps_point0 = mat_temps_parcours; 
    temps_point0 = temps_point0(index0,:)+recordedData.offsetTime{1,source_index}-temps_decalage;

    % {
    figure (1), 
    % on d�termine l'�chelle de temps des donn�es RF trait�es
    t0_donnees = recordedData.offsetTime{1,source_index}-temps_decalage; 
    taille_donnees_temp = size(recordedData.rfSignals{1,source_index},1); 
    vect_temps = t0_donnees:1/paramField.fs:(t0_donnees+(taille_donnees_temp-1)/paramField.fs);
    % on dessine ces signaux RF avec l'�chelle de temps d'aller retour
    imagesc(paramEmission.position_elements*10^3, vect_temps*10^6, recordedData.rfSignals{1,source_index}); 
    hold on
    % on dessine les temps utilis�s pour le diffuseur
    plot(paramEmission.position_elements(paramEmission.recepteurs_actifs)*10^3,temps_point0(paramEmission.recepteurs_actifs)*10^6,'+r'); %
    title(['Points utilis�s pour reconstruire l''amplitude du diffuseur']);%
    xlabel('Axe lat�ral de la sonde [mm]')
    ylabel('Temps [\mus]'); 
    colorbar
    
    figure (2), % pour voir le point suivi dans l'image basse r�solution reconstruite
    im_env= abs(hilbert( rfBeamformed(:,:,tir))); 
    im_env_norm_log = 20*log10(im_env/max(max(im_env)));  
    imagesc(x_grille*10^3,z_grille*10^3,im_env_norm_log), hold on, plot(x_grille*10^3, [zeros(1,num0-1) z_grille(numm0)*10^3 zeros(1,length(x_grille)-num0)],'xw', 'MarkerSize',24); 
    colormap(gray), caxis([-60 0]); 
    title([ 'Image basse r�solution, obtenue avec l''�mission  ' num2str(tir) ]); % / point utilise : x = ' num2str(imageGridX(num0)*10^3) ', z = ' num2str(imageGridZ(numm0)*10^3)]); 

    colormap(gray), caxis([-60 0]); 
    c = colorbar;
    c.Label.String = 'dB';
    c.Label.Position=[0.5 -60.8 0];
    c.Label.Rotation = 0;
    xlabel('Axe lat�ral de la sonde [mm]'); ylabel('Profondeur [mm]'); 

    %}
   % on stocke les temps utilis�s pour interpoler l'amplitude du point
   % suivi
   stockage_temps(:,tir) = temps_point0(paramEmission.recepteurs_actifs)';
   
   % on incr�mente le nb de tirs
   tir = tir +1; 
    end
delete(f);

end

