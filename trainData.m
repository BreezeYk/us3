    % %   main program implementing a simultaneous emission scheme in STA
    % utilizing diverging waves (code lengths can be tuned).
    % % the code generates rf signals 
    % % the offset due to the transducer convolution is taken into account 
    % % signal compression is performed before beamforming so that there is
    % no code length issues during dynamic focusing 
    
    %%
    clear all 
    close all 
    field_init() ;
    
    %%
    

    
    adresse_resultats = 'C:\Users\guinebretiere\Desktop\test'; 
    
    % ajustment / tuning / setting of the emission parameters
    
    paramEm.perfImpRep = 1 ;   % if perfImpRep = 1 : reponse impulsionnelle "parfaite" if 0 : reponse de la sonde Ulaop
    paramEm.exType = 1 ; % if exType = 1 : reponse avec exc codee , else dirac 
    paramEm.Na = 64 ;
    paramEm.nActiveArr = [ 1 : paramEm.Na ] ;
    paramEm.plotRep = 0 ;
    
    
    
    % % field parameters 
    
    paramField.mettre_chirp = 0 ; 
    paramField.sources_actives = paramEm.nActiveArr ; 
%     paramField.N_STA = length(paramField.sources_actives); % nb d'émissions / de tirs en ouverture synthétique
    paramField.N_STA = 1 ;
    paramField.fs = 100e6 ; % frequence d'echantillonnage
%     paramField.fs = 20832000 ; % frequence d'echantillonnage
    paramField.Ts = 1/paramField.fs;

    paramField.f0 = 8.5e6; 

    if paramEm.perfImpRep == 1
        paramField.f0 = 6e6;    % centre de la sonde
    end
    
    paramField.c = 1540;                                       % vitesse du son
    paramField.lambda = paramField.c/paramField.f0;            % longueur de l'onde d'emission
    paramField.hauteur_element = 6/1000;                       % Hauteur d'un element [m]
    paramField.kerf = 0.03/1000;                             % Distance entre les elements [m]
    paramField.largeur_element = 0.215/1000;                 % Largeur d'un element [m]
    paramField.pitch = paramField.kerf + paramField.largeur_element;  % pitch
    paramField.focus = [0 0 18]/1000;                        % Point de focus (x,y,z) [m]
    paramField.N_elements_emission = 64 ;                    % Nombre d'elements dans la barrette d emission
    paramField.N_elements_reception=paramField.N_elements_emission; % Nombre d'elements dans la barrette de reception
    paramField.retard = zeros(1,paramField.N_elements_reception); 
    paramField.adresse_resultats = adresse_resultats; 

    % % excitation signal 
    
    % generation of orthogonal Gold codes 
    
    paramEm.n = length( [1 1 0 0 1 1] ) ;
    paramEm.n5 = length( [1 0 0 1 0] ) ;
    paramEm.n9 = length( [1 0 0 0 0 1 0 0 0]  ) ;
    paramEm.n2 = length( [ 1 0 0 0 0 0 0 0 0 1 0 ] ) ;
    paramEm.L = 2^paramEm.n  ;
    paramEm.L5 = 2^paramEm.n5  ;
    paramEm.L2 = 2^paramEm.n2  ;
    paramEm.L9 = 2^paramEm.n9  ;
    
    if exist('g') == 0
        
        
        [g,G,l,OG5] = gold(  [1 0 0 1 0] , [1 1 1 1 0] , 1 ) ; % code of length L and 2^n sequences
        [g,G,l,OG] = gold(  [1 0 0 0 0 1] , [1 1 0 0 1 1] , 1 ) ; % code of length L and 2^n sequences
%         [OGt , l] = lkasami(  [1 0 0 0 0 1] , [1 1 0 0 1 1]  ) ; % code of length L and 2^n sequences
        [g9,G,l,OG9] = gold(  [1 0 0 0 0 1 0 0 0] , [1 0 0 1 0 1 1 0 0] , 1 ) ;
        [g2,G,l,OG2] = gold(  [ 1 0 0 0 0 0 0 0 0 1 0 ] , [ 1 0 0 1 0 0 1 0 0 1 0 ] , 1 ) ; % code of length L and 2^n sequences
        OG = 2*OG - 1 ;
        OG2 = 2*OG2 - 1 ;
        OG5 = 2*OG5 - 1 ;
        OG9 = 2*OG9 - 1 ;
        
    
    
    % BPSK with the Gold orthogonal codes 
    
    paramEm.excitation = sin( 2*pi*paramField.f0*(0:1/paramField.fs: 5/paramField.f0) ) ;  % carrier signal 
    paramEm.NbCycles = 1 ;  % number of periods per character of codes
    paramEm.nS = length( paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0)) ) ;  
    
    exArray = BPSK( OG , paramEm  , paramField , paramEm.L );
  %  exArray5 = BPSK( OG5 , paramEm  , paramField );
    exArray2 = BPSK( OG2 , paramEm  , paramField , paramEm.L2 );
    exArray9 = BPSK( OG9 , paramEm  , paramField , paramEm.L9  );
    
    clear OG OG2 g2 l g9 OG9;
    
    end
%     
%      ac1 = corr(exArray(1,:) , exArray(1,:))
%      cc1 = corr(exArray(1,:) , exArray(2,:))
%      
%      ac2 = corr(exArray2(1,:) , exArray2(1,:))
%      cc2 = corr(exArray2(1,:) , exArray2(2,:))
%     
    
    
% ex = direcGene( paramField.fs ) ;
    
    paramField.excitation = paramEm.excitation;
    paramModeleDirect.excitation = paramEm.excitation; 

    % % impulse response of the probe 
    
    load('LA523ERIat61p6MHz.mat');  % we load the Ula-op impulse response 
    
    if paramField.fs~=61.6e6        
    
        if size(impulse_response,1)==1
            
            impulse_response = impulse_response'; 
            
        end
        
            [impulse_response,~] = resample_fonction_MATRIX(impulse_response,61.6e6,paramField.fs); 
            impulse_response = impulse_response'; 
       
    end
    
    
    if paramEm.perfImpRep == 1 
        
        impulse_response_1 = paramEm.excitation ;
        impulse_response = impulse_response_1.*hanning(max(size(impulse_response_1)))' ; 
        
    end
    
    rep_impuls = impulse_response; 
%     clear impulse_response
    paramField.rep_impuls = rep_impuls ;
    
    
    if paramEm.plotRep == 1
        n = length(impulse_response);          % number of samples
        f = (0:n-1)*(fs/n);     % frequency range
        power = abs(fft(impulse_response)).^2/n;    % power of the DFT

        plot(f,power)
        xlabel('Frequency')
        ylabel('Power')
    end
   % % FIELD II
   
    field_init(0); 

    set_field('use_att', 0);              % ne pas utiliser de l'attenuation dans le calcul
    set_field('c',paramField.c);           % Initialisation de la vitesse
    set_field('fs',paramField.fs);         % Initialisation de la frequence d'echantillonage

    % création des ouvertures en émission et en réception 
    R=6 ;
%    for k = 1 : R
    pos = 3 ; 
    
        for p = 1:pos

        % % fantom => 2 scatterers
        
        z = 30:60 ;
        idxz = randi(length(z),1,1) ;
        x = -1:0.1:1 ;
        idxx = randi(length(x),1,1) ;
        
        pt_positions = [x(idxx)*1e-3 0 z(idxz)*1e-3 ]
        amp = [1]; 
       

        

        paramEm.exi = 6 ;

        if paramEm.exi == 6

            x = exArray ;
            
            mat_emission = zeros( paramField.N_elements_emission , length( x(1,:) ) );
            mat_emission(1,:) = x(1,:) ;
            mat_emission( floor( paramEm.Na*(1/4) )  ,:) = x(2,:) ;
            mat_emission( floor( paramEm.Na/2 ) , : ) = x(3,:) ;
            mat_emission( floor( paramEm.Na*(3/4) ) , : ) = x(4,:) ;
            mat_emission( end , :) = x(5,:) ;
            K = 64 ;
            for k = 1 : K
                
                mat_emissionT = zeros(K , paramField.N_elements_emission , length( x(1,:) ) ) ;
                mat_emissionT(k,1,:) = x(k,:) ;
                
            end
            mat_em = zeros( paramField.N_elements_emission , length( x(1,:) ) ) ;
            mat_em(1,:) = x(1,:) ;
            % stockage de cette matrice
            paramField.signalExcitation = exArray ;
            
        end

        
        for k = 1:K
        
            [emit_aperture receive_aperture] = aparture( paramField ) ;
            
            if k == 1 

                paramField.signalExcitationSimultane{1} = mat_emission ;

            else

             %  paramField.signalExcitationSimultane{1} = squeeze( mat_emissionT(k,:,:) ) ;
               paramField.signalExcitationSimultane{1} = mat_em ;
                
            end


            if paramEm.exType == 1 

                    ele_waveform(emit_aperture,(1:paramField.N_elements_emission)',paramField.signalExcitationSimultane{1}); 

            end

            if paramEm.exType == 0 

                duree_chirp = 5e-6; % durée en seconde
                largeur_bande_passante = paramField.f0; 
                [ex, ~, ~] = genChirp(duree_chirp, largeur_bande_passante, paramField);

        %         mat = zeros(16,201) ;
                mat = zeros( 16 , length( ex ) ) ;
                mat(1,:) = ex ;
                ele_waveform(emit_aperture,(1:paramField.N_elements_emission)', mat ); 

            end


            % calcul des reflexions. Cette fonction permet d'emettre les
            % excitations et recevoir les echos
            
            if k == 1
                [recordedData.rfSignals{1}, recordedData.offsetTime{1}] = calc_scat_multi( emit_aperture , receive_aperture, pt_positions, amp ); 
            end
            [ y(k,:,:) b ] = calc_scat_multi( emit_aperture , receive_aperture, pt_positions, amp ); 
            
            paramField.signalExcitationSimultane{1} ;
            
            clear emit_aperture
            clear receive_aperture
         %   clear paramField.signalExcitationSimultane{1}
        end
        
        for i = 1 : paramEm.Na    
         
            result( :,i ) = conv( conv( x(i,:) , paramField.rep_impuls ) , paramField.rep_impuls ); 
          
        end
        
        a = result( :,1 ) ;
        b = a(end:-1:1);
        yc = filter(b , 1 , squeeze( y(1,:,1) ) );
        
        
     %   figure() ; plot( y(k,1) )      
     a = y(1,:,1) ; 
     b = y(1:end,:,1) ; 
    save( ['train/mixed' num2str(p) '.mat' ] , 'a' ) ;
    ['train/mixed' num2str(p) '.mat' ]
    save( ['train/target' num2str(p) '.mat' ] , 'b' ) ;
    
    Lo = length( squeeze( y(1,:,1) ) )
    
    K = 5 ;
    X = zeros(K ,  Lo -  length(x(1,:)) , Lo ) ;
       
        for k = 1:K
            
            for i = 1 : Lo  - length(x(1,:))
                
               if i == 1
                   
                    X( k , i , i : i + length( x(1,:) ) -1  ) = x(k,:);
                    
               else
                    X( k , i , i : i + length( x(1,:) ) - 1  ) = x(k,:);
               end
            end
        end
        
        
    xk = squeeze( X(1,:,:) )' ;  
    yq = squeeze( y(1,:,1) ) ;
    h_hat = (inv( xk'*xk )*(xk'))*yq';
               
 %   y_hat = conv( x(1,:) , h_hat ) ;
    
    y_hat = squeeze( X(1,:,:) )' * h_hat  ;
   
 %  y_hat = conv( result(:,1) , h_hat )  ;
     sigrf = recordedData.rfSignals{1} ;
    [signauxRF1.rfSignals{1} , signauxRF2.rfSignals{1}] = compression_chirp_filtre_adapte( sigrf(:,1) , result(:,1) ); %excitation); 
    [signauxRF1.rfSignals{i} , yy] = compression_chirp_filtre_adapte( y_hat , result(:,1) ); %excitation); 
    
    figure() ; subplot(221) ; plot( squeeze( y(1,:,1) ) ) ; ylim([-0.1 0.1]*2e-23) ; subplot(222) ; plot( squeeze( y(2,:,1) ) ) ; ylim([-0.1 0.1]*2e-23)
    subplot(224) ; plot( y_hat ); ylim([-0.2 0.2]*1e-23)  ;
    
    r = signauxRF2.rfSignals{1} ; 
    figure() ; plot( 20*log10( abs( hilbert( r )))  ) ;
    figure() ; plot( 20*log10( abs( hilbert( yy )))  ) ;
    
    % pt_positions et amp viennent des données du fantôme
       clear emit_aperture
       clear receive_aperture
       clear y b ;
       
        end

        
    %% reconstruction DAS

% on récupère les paramètres d'émission
paramDAS = paramField; 
paramDAS.impulse_response = paramField.rep_impuls; 

% on calcule le resultat de la double convolution (excitation/transducteur)
clear result ; clear temps_result ; clear tmp ; clear max_pic ; clear tmp2 result_chirp_env2 max_pic2 temps_decalage3 ;
if paramEm.exType == 1
    for i = 1 : paramEm.Na    
        
        
        if paramEm.exi == 5
            
            result( :,i ) = conv( conv( exArray(i,:) , paramDAS.impulse_response ) , paramDAS.impulse_response ); 
            
        end   
        
        if paramEm.exi == 6
            
            result( :,i ) = conv( conv( exArray(i,:) , paramDAS.impulse_response ) , paramDAS.impulse_response ); 
            
        end    
        if paramEm.exi == 9
            
            result( :,i ) = conv( conv( exArray9(i,:) , paramDAS.impulse_response ) , paramDAS.impulse_response );
        end
        if paramEm.exi == 11
        
            result( :,i ) = conv( conv( exArray2(i,:) , paramDAS.impulse_response ) , paramDAS.impulse_response );
            
        end
    %     if size(result,1)>1
    %         result = result'; 
    %     end
        % on trouve le temps correspondant au max de l'enveloppe du signal
        temps_result(i , :) = (0:length(result(:,i))-1)/paramField.fs; 
        tmp(i , :)  = abs( hilbert( result(:,i) ) ); 
        max_pic(i)  = find( tmp( i , : ) == max( tmp( i , : ) ) ); 

    end
else 
    for i = 1 : 16    

        result( :,i ) = conv( conv( ex , paramDAS.impulse_response ) , paramDAS.impulse_response ); 
      
    %     if size(result,1)>1
    %         result = result'; 
    %     end
        % on trouve le temps correspondant au max de l'enveloppe du signal
        temps_result(i , :) = (0:length(result(:,i))-1)/paramField.fs; 
        tmp(i , :)  = abs( hilbert( result(:,i) ) ); 
        max_pic(i)  = find( tmp( i , : ) == max( tmp( i , : ) ) ); 

    end
    
end
% parametres de l'image à reconstruire
paramImage.prof_min_milieu = 40e-3 ; % en mètre
paramImage.prof_max_milieu = 60e-3 ; % en mètre
paramImage.resolution_axiale = paramField.c/(2*paramField.fs); % en mètre
paramImage.x_min = -5e-3 ; % en mètre
paramImage.x_max = 5e-3; % en mètre
paramImage.resolution_laterale = paramDAS.position_elements(1,2)-paramDAS.position_elements(1,1); % en mètre
paramImage.position_a_suivre = pt_positions(1,:); % mettre les coordonnées 
% d'un diffuseur ou d'un point quelconque (obligatoire), cela permet de 
% suivre ce point pendant la reconstruction
paramImage.diffuseur_reso = ''; % mettre les coordonnées d'un point ou ''. 
% Cette option permet de reconstruire l'image en étant centré sur ce point 
% (utile pour le calcul de la résolution)


% on comprime les signaux RF par leurs filtres adaptés respectifs
wien = 0 ;
if 1
    for i = 1 : paramEm.Na
        
        if wien == 0 
            
            tmp2(i,:) = compression_chirp_filtre_adapte( result(:,i) , result(:,i) );
            
        else
        %    tmp2(i,:) = compression_chirp_filtre_wiener3( result(:,i) , result(:,i) , 25);
        end
        result_chirp_env2(i,:) = abs(hilbert(tmp2(i,:))); 
        [~, max_pic2(i)] = max(result_chirp_env2(i,70:400));
%         [~, max_pic2(i)] = max(result_chirp_env2(i,:));
        temps_decalage3(i) = temps_result( max_pic2(i) + 70 ); % temps de décalage de l'enveloppe
        
%             for source_index = paramDAS.sources_actives
%                 disp(['Compression des chirps ' num2str(source_index) ' sur ' num2str(paramField.N_STA)]); 

                % on comprime les signaux RF avec adapté
                sigrf = recordedData.rfSignals{1} ; 
                offset = recordedData.offsetTime{1} ;
                if wien == 0 
                    [signauxRF1.rfSignals{i} , signauxRF2.rfSignals{i}] = compression_chirp_filtre_adapte( sigrf(:,i) , result(:,5) ); %excitation); 
                else
             %       [signauxRF1.rfSignals{i} , signauxRF2.rfSignals{i}] = compression_chirp_filtre_wiener3( sigrf(:,i) , result(:,3) , 25 ); %excitation); 
                end
                signauxRF2.offsetTime{i} = offset - (length( result(i,:) )-1)/paramField.fs;  
                signauxRF1.offsetTime{i} = offset ; 
%             end
            
            adresse2 = [adresse_resultats '/adapte/']; 
        
        temps_decalage_final(i) = temps_decalage3(i) ; 
        if exist(adresse2)==0
            mkdir(adresse2); 
        end
%     else % si pas de compression
%         temps_decalage_final = temps_result(max_pic); %length(result)/2/paramField.fs; 
%         signauxRF2 = recordedData;
%         adresse2 = [adresse_resultats '/chip/']; 
%         if exist(adresse2)==0
%             mkdir(adresse2); 
%         end
end
end

    