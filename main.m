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
    
    paramEm.perfImpRep = 0 ;   % if perfImpRep = 1 : reponse impulsionnelle "parfaite" if 0 : reponse de la sonde Ulaop
    paramEm.exType = 1 ; % if exType = 1 : reponse avec exc codee , else dirac 
    paramEm.Na = 64 ;
    paramEm.nActiveArr = [ 1 : paramEm.Na ] ;
    paramEm.plotRep = 1 ;
    paramEm.K=3; % number of different emission codes
    
    
    
    % % field parameters 
    
    paramField.mettre_chirp = 0 ; 
    paramField.sources_actives = paramEm.nActiveArr ; 
%     paramField.N_STA = length(paramField.sources_actives); % nb d'émissions / de tirs en ouverture synthétique
    paramField.N_STA = 1 ;
    paramField.fs = 200e6 ; % frequence d'echantillonnage
%     paramField.fs = 20832000 ; % frequence d'echantillonnage
    paramField.Ts = 1/paramField.fs;

    paramField.f0 = 8.5e6; 

    if paramEm.perfImpRep == 1
        paramField.f0 = 6e6;    % centre de la sonde
    end
    
    
    
    
    paramField.c = 1540;                                       % vitesse du son
    paramField.lambda = paramField.c/paramField.f0;            % longueur de l'onde d'emission
    paramField.hauteur_element = 6/1000;
%     paramField.largeur_element = 1/1000;% Hauteur d'un element [m]
    paramField.kerf = 0.03/1000;   
%     paramField.kerf = paramField.largeur_element/5;   % Distance entre les elements [m]
    paramField.largeur_element = 0.215/1000;                 % Largeur d'un element [m]
    paramField.pitch = paramField.kerf + paramField.largeur_element;  % pitch
    paramField.focus = [0 0 18]/1000;                        % Point de focus (x,y,z) [m]
    paramField.N_elements_emission = paramEm.Na ;                    % Nombre d'elements dans la barrette d emission
    paramField.N_elements_reception=paramField.N_elements_emission; % Nombre d'elements dans la barrette de reception
    paramField.retard = zeros(1,paramField.N_elements_reception); 
    paramField.adresse_resultats = adresse_resultats; 

    % % excitation signal 
    
    
    paramEm.exi = 6 ;
   % ex = direcGene( paramField.fs , size(x6,2) ) ;
   
   

    % % impulse response of the probe 
    
    load('LA523ERIat61p6MHz.mat');  % we load the Ula-op impulse response 
    
    if paramField.fs~=61.6e6        
    
        if size(impulse_response,1)==1
            
            impulse_response = impulse_response'; 
            
        end
        
            [impulse_response,~] = resample_fonction_MATRIX(impulse_response,61.6e6,paramField.fs); 
            impulse_response = impulse_response'; 
       
    end
    
  
    
    rep_impuls = impulse_response; 
%     clear impulse_response
    paramField.rep_impuls = rep_impuls ;
    
    
    if paramEm.plotRep == 1
        n = length(impulse_response);          % number of samples
        f = (0:n-1)*(paramField.fs/n);     % frequency range
        power = abs(fft(impulse_response)).^2/n;    % power of the DFT
        figure; 
        plot(f,power)
        xlabel('Frequency')
        ylabel('Power')
    end
    
    % generation of orthogonal Gold codes 
    
    paramEm.n = length( [1 1 0 0 1 1] ) ;
    paramEm.n5 = length( [1 0 0 1 0] ) ;
    paramEm.n9 = length( [1 0 0 0 0 1 0 0 0]  ) ;
    paramEm.n7 = 7 ;
    
    paramEm.n10 = 10 ;
    paramEm.n2 = length( [ 1 0 0 0 0 0 0 0 0 1 0 ] ) ;
    paramEm.L = 2^paramEm.n  ;
    paramEm.L5 = 2^paramEm.n5  ;
    paramEm.L2 = 2^paramEm.n2  ;
    paramEm.L9 = 2^paramEm.n9  ;
    paramEm.L10 = 2^paramEm.n10  ;
    
    if exist('g') == 0
        
        
        [g,G,l,OG5] = gold(  [1 0 0 1 0] , [1 1 1 1 0] , 1 ) ; % code of length L and 2^n sequences
        [g,G,l,OG] = gold(  [1 0 0 0 0 1] , [1 1 0 0 1 1] , 1 ) ; % code of length L and 2^n sequences
%         [OGt , l] = lkasami(  [1 0 0 0 0 1] , [1 1 0 0 1 1]  ) ; % code of length L and 2^n sequences
        [g9,G,l,OG9] = gold(  [1 0 0 0 0 1 0 0 0] , [1 0 0 1 0 1 1 0 0] , 1 ) ;
        [g7,G,l,OG7] = gold(  [1 0 0 0 1 0 0] , [1 0 0 0 1 1 1] , 1 ) ;
        [g10,G,l,OG10] = gold(  [1 0 0 0 0 0 0 1 0 0] , [1 0 1 0 0 0 0 1 1 0] , 1 ) ;
        [g2,G,l,OG2] = gold(  [ 1 0 0 0 0 0 0 0 0 1 0 ] , [ 1 0 0 1 0 0 1 0 0 1 0 ] , 1 ) ; % code of length L and 2^n sequences
        OG = 2*OG - 1 ;
        OG2 = 2*OG2 - 1 ;
        OG5 = 2*OG5 - 1 ;
        OG9 = 2*OG9 - 1 ;
        OG10 = 2*OG10 - 1 ;
        OG7 = 2*OG7 - 1 ;
        
    
   
    % BPSK with the Gold orthogonal codes 
    
    paramEm.excitation = sin( 2*pi*paramField.f0*(0:1/paramField.fs: 5/paramField.f0) ) ;     % carrier signal 
    paramEm.NbCycles = 2 ;  % number of periods per character of codes
   % paramEm.nS = length( paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0)) ) ;  
    paramEm.nS = length( sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) ) ) ;
    
    binary = 1 ;
    sparse = 0 ;
    if sparse == 1 
        win = 0 ;
        %x6 = BPSK1( OG , paramEm  , paramField , paramEm.L );
        x6 = BPSK1( OG , paramEm  , paramField  , 1 , win , binary );
        x7 = BPSK1( OG7 , paramEm  , paramField  , 4 , win, binary);
    %     exArray5 = BPSK( OG5 , paramEm  , paramField , paramEm.L5 );
        x11 = BPSK1( OG9 , paramEm  , paramField , 4 , win, binary);
        x9 = BPSK1( OG9 , paramEm  , paramField  , 4 , win, binary);
        x10 = BPSK1( OG7 , paramEm  , paramField , 4 , win, binary);
        
       
    else
        win = 0 ;
        %x6 = BPSK1( OG , paramEm  , paramField , paramEm.L );
        x6 = BPSK( OG , paramEm  , paramField , paramEm.L , win  );
        x7 = BPSK( OG7 , paramEm  , paramField , paramEm.L , win);
    %     exArray5 = BPSK( OG5 , paramEm  , paramField , paramEm.L5 );
        x11 = BPSK( OG2 , paramEm  , paramField ,  paramEm.L2 , win);
        x9 = BPSK( OG9 , paramEm  , paramField ,  paramEm.L9 , win);
        x10 = BPSK( OG10 , paramEm  , paramField ,  paramEm.L10 , win);
    
    end
%     clear OG OG2 g2 l g9 OG9;
    
%     
%      ac1 = corr(x6(1,:) , x6(1,:))
%      cc1 = corr(x6(1,:) , x6(2,:))
%      
%      ac2 = corr(x11(1,:) , x11(1,:))
%      cc2 = corr(x11(1,:) , x11(2,:))
%     
    end
    
        binary = 1 ;
        g = normrnd(0, 1 , 20 , 512) ;
        
        a = 2*(g>0) - 1 ;
        
        win = 0 ;
        
     %   xx = eval(['x' num2str(paramEm.exi)]) ; 
      %  xx = BPSK1(a , paramEm  , paramField  , 2 , win , 0) ;
        xx = BPSK1(OG9 , paramEm  , paramField  , 2 , win , 0) ;
       % xx1 = BPSK1(a , paramEm  , paramField  , 2 , win , 1) ;
        xx1 = BPSK1(OG9 , paramEm  , paramField  , 2 , win , 1) ;
    
        % % Waveform chip 
        
        chip = sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) ) ; 
        chipp = conv( conv( chip ,  paramField.rep_impuls  ) , paramField.rep_impuls ) ;
        
      %  xx = x6 ;
     ex = direcGene( paramField.fs , size(xx,2) ) ;
    
    paramField.excitation = paramEm.excitation;
    paramModeleDirect.excitation = paramEm.excitation; 
    
    if paramEm.perfImpRep == 1 
        
        impulse_response_1 = paramEm.excitation ;
        impulse_response = impulse_response_1.*hanning(max(size(impulse_response_1)))' ; 
        rep_impuls = impulse_response; 
%     clear impulse_response
        paramField.rep_impuls = rep_impuls ;
        
    end
    
    clear Y Yc C h_hat x y_hat target obs yt train h_target
    
    nPos = 1
    tic
    for pos = 1 : nPos
    
        pos
    
    
   % % FIELD II
   
    field_init(0); 
    
    set_field('use_att', 0);              % ne pas utiliser de l'attenuation dans le calcul
    set_field('c',paramField.c);           % Initialisation de la vitesse
    set_field('fs',paramField.fs);         % Initialisation de la frequence d'echantillonage

    % création des ouvertures en émission et en réception 
%     emit_aperture = xdc_linear_array (paramField.N_elements_emission,...
%                                       paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
%                                       1, 1,paramField.focus);    % le focus qui est mis ici est temporaire et ne va pas servir a la focalisation, 
%                                                                  % car apres en imposant les delays (signaux d'excitation) on supprime l'effet de ce focus
% 
% 
%                                                                  emit_apertureT = xdc_linear_array (paramField.N_elements_emission,...
%                                       paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
%                                       1, 1,paramField.focus); 
%                                                                  
    receive_aperture = xdc_linear_array (paramField.N_elements_reception,...
                                      paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
                                      1, 1,paramField.focus);    % le focus qui est mis ici est temporaire et ne va pas servir a la focalisation, 
                                                                 % car apres en imposant les delays (signaux d'excitation) on supprime l'effet de ce focus

    
    % % 

    % % fantom => 2 scatterers
    
    
    xpos = [-2 : 2]*1e-3;
    zpos = [30 : 60]*1e-3;
    
    idxx = randperm( length(xpos) ); 
    idxz = randperm( length(zpos) );
        
        scat = 1 ;

        if scat == 2 
            amp = [1;1;1];
          %  pt_positions = [0 0 25e-3 ; -1e-3 0 30e-3 ;  1e-3 0 45e-3];
            pt_positions = [0 0 25e-3 ; 0 0 30e-3 ;  0 0 35e-3];
        elseif scat == 1
            amp = [1;1];
            pt_positions = [0 0 20e-3 ; -1e-3 0 22e-3 ];
        else
            pt_positions = [0 0 25e-3 ]; 
            amp = [1];
           % pt_positions = [xpos( idxx(1) ) 0 zpos( idxz(1) ) ] 
           % amp = [1];
        end
         
        arrayEm = [1 2 33 34] ;
       paramEm.K = length( arrayEm ) ;
    for k = 1 : 2*(paramEm.K)+1
        
            emit_aperture(k) = xdc_linear_array (paramField.N_elements_emission,...
                                      paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
                                      1, 1,paramField.focus); 
        
            paramField.probe_specs = xdc_get(emit_aperture(k) , 'rect');
%             paramField.probe_specs = xdc_get(emit_apertureT , 'rect');% on recupere tous les parameters geometriques de la probe
            paramField.position_elements=paramField.probe_specs(8,:); % on recupere les coordonees des centres des elements de la probe
           
            xdc_impulse(emit_aperture(k), rep_impuls);      % en emission
            xdc_impulse(receive_aperture, rep_impuls);   % en reception

            % On supprime le focus en emission
            xdc_center_focus(emit_aperture(k),[0,0,0]);
            xdc_center_focus(receive_aperture,[0,0,0]);

            % On force les delays d'emission a 0 dans un premier temps
            xdc_focus_times(emit_aperture(k), 0, zeros(1,paramField.N_elements_emission));        % en emission 
            xdc_focus_times(receive_aperture, 0, zeros(1,paramField.N_elements_reception));    % en reception

            % On supprime le focus en emission
            xdc_center_focus(emit_aperture(k),[0,0,0]);

            xdc_apodization( emit_aperture(k),0, ones(1,paramField.N_elements_emission) );
            xdc_apodization( receive_aperture,0, ones(1,paramField.N_elements_reception) );
    
            end 
            
    chip = paramField.excitation( 1:floor(2.5*paramField.fs/paramField.f0) ) ;
    
    paramEm.h = tukeywin( size(chip,2) , 0.8 ) ;
    chip = ( chip .* paramEm.h' ) ;
     

    if paramEm.exi == 0
       
        mat_emission = zeros( paramField.N_elements_emission , length( chip ) );
        mat_emission(32,:) = chip ;
        paramField.signalExcitation = chip ;
        
    end
     
    %   arrayEm = 1:10:paramEm.Na ;
     %  arrayEm = [1 15 33 47] ;
       
       paramEm.K = length( arrayEm ) ;
       paramEm.K
       mat_emission = zeros( 2*length(paramEm.K)+1 , paramField.N_elements_emission , length( xx(1,:) ) );
       
       r=1;
       
      for j = arrayEm
          
         mat_emission(1,j, : ) = xx(r,:) ;
         mat_emission(r+1,j,:) = [chip zeros( length( xx(1,:) ) - length(chip) , 1 )' ] ;
           %  mat_emission(r+1,j,:) = ex ;
         mat_emission(r+1+paramEm.K ,j,:) = xx(r,:) ;
         
         r = r+1;
         
      end

         % stockage de cette matrice
        paramField.signalExcitation = xx ;
     
    
%     emit_apertureT = emit_aperture ;
    
    if paramEm.exType == 1 
    
        for k = 1 : 2*(paramEm.K)+1
            
          paramField.signalExcitationSimultane{k} = squeeze( mat_emission(k,:,:) );
          ele_waveform(emit_aperture(k),(1:paramField.N_elements_emission)',paramField.signalExcitationSimultane{k});
%           ele_waveform(emit_apertureT ,(1:paramField.N_elements_emission)', mat_emissionTarget );

        end
        
    end
    
    if paramEm.exType == 0 
        
        duree_chirp = 5e-6; % durée en seconde
        largeur_bande_passante = paramField.f0; 
        [ex, ~, ~] = genChirp(duree_chirp, largeur_bande_passante, paramField);
        
%         mat = zeros(16,201) ;
        mat = zeros( paramEm.Na , length( ex ) ) ;
        mat(1,:) = ex ;
        ele_waveform(emit_aperture,(1:paramField.N_elements_emission)', mat ); 
    
    end
  

    % calcul des reflexions. Cette fonction permet d'emettre les
    % excitations et recevoir les echos
    
    for k = 1:2*(paramEm.K)+1
    
        
        if k ==1
            
            [recordedData.rfSignals{1} , recordedData.offsetTime{1}] = calc_scat_multi( emit_aperture(k) , receive_aperture, pt_positions, amp ); 
        l=1 ;
        elseif ( k < paramEm.K+2 && k > 1 )
            k
            [h_target(l,:,:) , os] = calc_scat_multi( emit_aperture(k) , receive_aperture, pt_positions, amp ); 
            l = l + 1 ;
            
        l1 = 1 ;       
        elseif k >= paramEm.K+2
            k
            [yt(l1,:,:) b] = calc_scat_multi( emit_aperture(k) , receive_aperture, pt_positions, amp ); 
            l1 = l1 + 1 ;
        end
    end
    
        Yc = recordedData.rfSignals{1} ;
    
        X = zeros( paramEm.K , size( Yc , 1 ) , size( Yc , 1 ) - size( xx , 2 ) + 1  ) ;
        M =     size( Yc , 1 ) - size( xx , 2 ) + 1 ;
        x = xx1 ;
     %   paramEm.K=6;
    
 %for nM = 15:10:60   
 %   for nM = 55:55  
    clear h_targetN
     nM = 5 ;
    for i = 1 : paramEm.Na
        
       Y( : , i ) = add_noise( Yc( : , i ) , nM , paramField , '' ) ;
       for k = 1 : paramEm.K
            h_targetN( : , i ) = add_noise( squeeze( h_target(k, : , i ) ) , nM , paramField , '' ) ;
       end
       
    end
    
    paramEm.Kt = length(arrayEm) ;
 %   paramEm.Kt = 3 ;
 %   for k = 1:size(arrayEm,2)
    for k = 1:paramEm.Kt
        
        for i = 1 : M 
            
            X( k , i : i + size( x , 2 ) - 1 , i ) = x(k,:) ;

        end
        
    end
    
    Xt=[];
    for k = 1:paramEm.Kt
        Xt = [Xt squeeze(X( k , : , : )) ];
    end
    
    size(Xt)
    
    if ( size(Xt,1) < size(Xt,2) )
       disp('fatal error') 
    end
    
    sep = 1 ;
    
    if sep == 1
    
    LL = 1 ;
    % L = 1;
    clear h_hatt
    
    
    for q = 1:LL
        
        q
        er = inv( mtimes( Xt' , Xt ) ) ;
        for k = 1 : paramEm.K
            for m = 1 : M    
                esErr(q,m,k) = er( k*(M)-M+m , k*(M)-M+m ) ;
            end
        end
        h_hatt(:,q) = mtimes( er ,  Xt') * Y(:,q) ;
      %  h_hatt(:,q) = ( mtimes( Xt' , Xt )  \  Xt' ) * Y(:,q) ;
        
    end
    
   % for k = 1:size(arrayEm,2)
    %for k = 1:paramEm.K
      %  k
      % for q = 1 : paramEm.Na
     %      q
        %    h_hat(k,q,:) =  mtimes( inv( mtimes( squeeze( X(k,:,:) )' , squeeze( X(k,:,:) ) ) ) , squeeze( X(k,:,:) )') * Yc(:,q)  ;
        %   y_hat(k,q,:) = conv( x(k,:) , squeeze( h_hat(k,q,:) ) ) ;
    %       y_hat(q,:) = conv( x(k,:) , squeeze( h_hatt(1:length(h_hatt)/3,q) ) ) ;
            
  %     end
        
  %  end
    
    end
    
    paramDAS = paramField; 
paramDAS.impulse_response = paramField.rep_impuls; 
clear result ; clear temps_result ; clear tmp ; clear max_pic ; clear tmp2 result_chirp_env2 max_pic2 temps_decalage3 ;
if paramEm.exType == 1 %result computation
    for i = 1 : paramEm.Na    
        
      for j = 1 : size(xx,1)
        
            result( :,i,j ) = conv( conv( xx(j,:) , paramDAS.impulse_response ) , paramDAS.impulse_response );
            
      end
    %     if size(result,1)>1
    %         result = result'; 
    %     end
        % on trouve le temps correspondant au max de l'enveloppe du signal
   %     temps_result(i , :) = (0:length(result(:,i))-1)/paramField.fs; 
   %     tmp(i , :)  = abs( hilbert( result(:,i) ) ); 
   %     max_pic(i)  = find( tmp( i , : ) == max( tmp( i , : ) ) ); 

    end
end
    
    
    
    isl = 0 ;
 
    
   if isl == 1
    
    s = squeeze( result(:,:,1) ) ;
    st = Yc ;
    
    N = size( s , 1 ) ; N1 = size( st , 1 ) ;
    p = N ;
    K = N + 2*p ;
    s = cat(1 , s , zeros(p,paramEm.Na) ); s = cat(1 , zeros(p,paramEm.Na) , s );
    st = cat(1 , st , zeros(H-N1,paramEm.Na) ); st = cat(1 , zeros(H-N1,paramEm.Na) , st );
    S = zeros( paramEm.Na , 2*K-1 , H ) ;
    St = zeros( paramEm.Na , 2*K-1 , H ) ;
    sf = flip( s ) ; sft = flip( st ) ;
       
       for k = 1: paramEm.Na
           for i = 1 : H 
                
                
               S( k , i : i + H - 1 , i ) = sf(:,k)  ;
               St( k , i : i + H - 1 , i ) = sft(:,k)  ;
                
           end
       end
   
    
    
   ar = ones(2*K-1,1) ;
   ar( N+p ) = 0 ;
   % F = diag( ar )*(1/(2*H-1)) ;
   F = diag( ar )*(1/(2*H-1)) ;
    
   B = zeros(H,H) ;
   B = ( squeeze(S(1,:,:))'*F*squeeze(S(1,:,:)) ) ;
   for i = 1 : 3
        
      B = B + squeeze(S(i,:,:))'*F*squeeze(S(i,:,:)) ;
        
   end
    
   q_isl = (( squeeze( s(:,1) )'*squeeze( s(:,1) ))*inv( squeeze(S(1,:,:))'*F*squeeze(S(1,:,:)) )*squeeze(s(:,1)))/((squeeze(s(:,1))')*inv( squeeze(S(1,:,:))'*F*squeeze(S(1,:,:)))* squeeze(s(:,1)) ) ;
  % q_isl = inv(B)* squeeze( s(:,1) )*inv(squeeze( s(:,1) )'*B*squeeze(s(:,1)))*N ;
    
   y_isl = squeeze(St(1,:,:))*q_isl ; 
  % y_isl = squeeze(St(1,:,:))*s(:,1) ; 
    
   figure(); plot( y_isl ) ;    
    
    end
 %  subplot(223) ; plot( squeeze( y_hat(q,:) ) ) ; ylim([-1 1]*1e-24)
    
  %  rmse = sum ( sqrt( abs( squeeze( y_hat(k,q,:) ) - squeeze( yt(1,:,q) ) ).^2 ) ) ;
  %   rmse = sum ( sqrt( abs( squeeze( Yc(:,q) ) - squeeze( yt(1,:,q) ) ).^2 ) ) ;
    
   % figure(); subplot(121) ; plot(20*log10( abs( hilbert( xcorr( squeeze( y_hat(k,q,:) ) , squeeze( yt(1,:,q) ) )) ) )) ; ...
   %     subplot(122) ; plot(20*log10( abs( hilbert ( xcorr( squeeze( Yc(:,q) ) , squeeze( yt(1,:,q) ) )))))
   
 TR = 0 ;   
 if TR == 1  
    
    for k = 1:paramEm.K
    
        meanyt = mean( squeeze( yt(k,:,q) ) ) ; 
        stdyt = std( squeeze( yt(k,:,q) ) ) ;
        target(k,:) = ( squeeze( yt(k,:,q) ) - meanyt ) / stdyt ;
    
    end
    
    meanyc = mean( Yc(:,q) ) ; 
    stdyc = std( Yc(:,q) ) ;
    obs = ( Yc(:,q) - meanyc ) / stdyc ;
    
    % pt_positions et amp viennent des données du fantôme
    clear emit_aperture
    clear receive_aperture
    
    for k = 1:paramEm.K+1
        if k == 1
            
            train(k,:) = obs ;
            
        else
            train(k,:) = target(k-1,:)';
            
        end
    end
    save( ['train' num2str(pos) '.mat'] , 'train' ) ;
    
    end
%     clearvars -except train
   if pos ~= nPos
       clear Y Yc C h_hat x y_hat target obs yt train
   end
 
%     
%     figure()
%     figure(); 
%     subplot(2,2,1) ; plot( cc1 ) ; ylim([-500 1000]) ; title('CC L = 64')
%     subplot(2,2,2) ;  plot( ac1 ) ; title('AC  L = 64 ')
%     subplot(2,2,3) ; plot( cc2 ) ; ylim([-1000 10000]) ; title('CC L = 2048')
%     subplot(2,2,4) ;  plot(ac2 ) ; title('AC L = 2048')


k = 1 ;
q = 1 ;
figure(); subplot(221) ; plot( squeeze( yt(1,:,q) ) ); title('y11'); ylim([-1 1]*1e-24) ; subplot(222) ; plot( Y(:,q) ) ; title('y1') ; ylim([-1 1]*1e-24) ; % ... % subplot(223) ; plot( y_isl ) ; 
subplot(223)
%plot( squeeze( yt(2,:,q) ) ); title('y12'); ylim([-1 1]*1e-24) ;

%  figure();  plot(20*log10( abs( hilbert ( xcorr( squeeze( Yc(:,q)  ) , squeeze( yt(1,:,q) ) ))))) ; ylim([-1000 -900]); %  ...


%      hold on ; legend('decoded' ,'undecoded')  ; 
 %  plot(20*log10( abs( hilbert( xcorr( squeeze( y_hat(k,q,:)  ) , squeeze( yt(1,:,q) ) )) ) )); 
    
%  figure();  plot(20*log10( abs( hilbert ( xcorr( squeeze( Yc(:,q)  ) , squeeze( result(:,1,1) ) ))))) ; ylim([-450 -350]);  ...
     %   hold on ; plot(20*log10( abs( hilbert ( xcorr( y_isl , result(:,1) ))))) ; 
  %      hold on ;  plot(20*log10( abs( hilbert( xcorr( squeeze( y_hat(k,q,:)  ) , result(:,1) )) ) )); legend('decoded' ,'undecoded')  ; 

toc;
toc-tic

% figure(); subplot(121); plot(h_target(:,1)) ; % xlim([0 300]) ;
%   subplot(122); plot( squeeze( h_hatt(1:length(h_hatt)/paramEm.Kt,1) ) ) ; title(['Poly char' num2str(paramEm.exi) '  ,'  num2str(toc)  ' sec ' ' nb sensors :' num2str(size(arrayEm,2)) ])

%for i = 1 : length(arrayEm)
for i = 1 : LL
    for k = 1 : paramEm.K
       if k == 1
           h_est(k,:,i) = squeeze( h_hatt( 1:length(h_hatt)/paramEm.K , i ) ) ;
       else
           h_est(k,:,i) = squeeze( h_hatt( (k-1)*length(h_hatt)/paramEm.K:(k)*length(h_hatt)/paramEm.K-1, i ) ) ;
       end
    end
end

figure();
l=1 ; 
for k = 1 : paramEm.K   
   subplot(paramEm.K,2,l); plot( squeeze( h_target(k,:,1)) ) ; title(['h target  ' num2str(k) '1']); xlim([0 length(h_hatt)/paramEm.K]) ;
   if k == 1
       subplot(paramEm.K,2,l+1); plot( squeeze( h_hatt( 1:length(h_hatt)/paramEm.K , 1 ) ) ) ;  title(['h est  ' num2str(k) '1 ' 'noise (dB) : ' num2str(nM) ]); xlim([0 length(h_hatt)/paramEm.K]) ;
   else
       subplot(paramEm.K,2,l+1); plot( squeeze( h_hatt( (k-1)*length(h_hatt)/paramEm.K:(k)*length(h_hatt)/paramEm.Kt, 1 ) ) ) ;  title(['h est  ' num2str(k) '1  ' 'noise (dB) : ' num2str(nM)]); xlim([0 length(h_hatt)/paramEm.K]) ; 
   end
   l=l+2;
end




for i = 1 : LL
    i
    for k = 1 : paramEm.K
        mf(i,:,k) = compression_chirp_filtre_adapte( squeeze( h_est(k,:,i) )' , chipp  ) ; 
    end
end
  %end
    end  
    %% reconstruction DAS

% on récupère les paramètres d'émission
paramDAS = paramField; 
paramDAS.impulse_response = paramField.rep_impuls; 

% on calcule le resultat de la double convolution (excitation/transducteur)
clear result ; clear temps_result ; clear tmp ; clear max_pic ; clear tmp2 result_chirp_env2 max_pic2 temps_decalage3 ;
if paramEm.exType == 1
    for i = 1 : paramEm.Na    
        
            result( :,i ) = conv( conv( xx(i,:) , paramDAS.impulse_response ) , paramDAS.impulse_response );
            
        
    %     if size(result,1)>1
    %         result = result'; 
    %     end
        % on trouve le temps correspondant au max de l'enveloppe du signal
     %   temps_result(i , :) = (0:length(result(:,i))-1)/paramField.fs; 
     %   tmp(i , :)  = abs( hilbert( result(:,i) ) ); 
     %   max_pic(i)  = find( tmp( i , : ) == max( tmp( i , : ) ) ); 

    end
else 
    for i = 1 : paramField.N_elements_emission    

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
wien = 1  ;
betha = 25 ;
q = 32 ;
k = 1 ;
if 1
    for i = 1 : paramEm.Na
        
        if wien == 0 
            tmp2(i,:) = compression_chirp_filtre_adapte( result(:,i) , result(:,i) );
        else
            tmp2(i,:) = compression_chirp_filtre_wiener3( result(:,i) , result(:,i) , betha);
        end
        result_chirp_env2(i,:) = abs(hilbert(tmp2(i,:))); 
        [~, max_pic2(i)] = max(result_chirp_env2(i,70:400));
%         [~, max_pic2(i)] = max(result_chirp_env2(i,:));
%        temps_decalage3(i) = temps_result( max_pic2(i) + 70 ); % temps de décalage de l'enveloppe
        
%             for source_index = paramDAS.sources_actives
%                 disp(['Compression des chirps ' num2str(source_index) ' sur ' num2str(paramField.N_STA)]); 

                % on comprime les signaux RF avec adapté
                sigrf = recordedData.rfSignals{1} ; 
                offset = recordedData.offsetTime{1} ;
                if wien == 0 
                    [signauxRF1.rfSignals{i} , signauxRF2.rfSignals{i}] = compression_chirp_filtre_adapte( sigrf(:,i) , result(:,k) ); %excitation); 
                 %   [matchedY1(i,:) , matchedY(i,:)] = compression_chirp_filtre_adapte( squeeze( y_hat(k,i,:) ) , result(:,k) ); %excitation); 
                else
                    [signauxRF1.rfSignals{i} , signauxRF2.rfSignals{i}] = compression_chirp_filtre_wiener3( sigrf(:,i) , result(:,k) , betha ); %excitation);
                 %   [matchedY1(i,:) ,matchedY(i,:)] = compression_chirp_filtre_wiener3( squeeze( y_hat(k,i,:) ) , result(:,k) , betha ); %excitation); 
                end
                signauxRF2.offsetTime{i} = offset - (length( result(i,:) )-1)/paramField.fs;  
                signauxRF1.offsetTime{i} = offset ; 
%             end
            
            adresse2 = [adresse_resultats '/adapte/']; 
        
%        temps_decalage_final(i) = temps_decalage3(i) ; 
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


figure(); subplot(223) ; plot( signauxRF2.rfSignals{q} ) ; subplot(224) ; plot( 20*log10( abs( hilbert( signauxRF2.rfSignals{q} ) ) ) )  ; subplot(221) ; plot(Yc(:,q)) ;
% figure(); subplot(223) ; plot( matchedY(q,:) ) ; subplot(224) ; plot( 20*log10( abs( hilbert( matchedY(q,:) ) ) ) )  ; subplot(221) ; plot(squeeze( y_hat(1,1,:) )) ;

% on rajoute des champs pour le DAS
% paramDAS.sources_emettrices = paramDAS.sources_actives(1:end); 
% paramDAS.sources_emettrices = 1; 
% paramDAS.recepteurs_actifs = [1:16]; 
% 
% % on reconstruit les images basse résolution par DAS
% [rfBeamformed, x_grille, z_grille, stockage_temps, num0,numm0, amplitude0] = reconstruire_image_STA_codes_recepteurs(signauxRF2,paramImage,paramDAS,paramDAS, temps_decalage_final(1)); 

% % on reconstruit l'image par somme de ces images basse résolution
% IMAGE = sum(rfBeamformed,3);
% 
% % on visualise le résultat 
% figure, imagesc(x_grille*10^3,z_grille*10^3 , IMAGE); title('Image reconstruite par ouverture synthetique')
% xlabel('Axe latéral de la sonde [mm]'); ylabel('Profondeur [mm]'); 
% 
% % Post-traitements (extraction de l'enveloppe + échelle log) : 
% IMAGE_env = abs(hilbert(IMAGE)); 
% IMAGE_env_norm = IMAGE_env/max(max(IMAGE_env)); 
% IMAGE_env_norm_log = 20*log10(IMAGE_env_norm); 
% 
% % Affichage de l'image mode B
% figure, 
% imagesc(x_grille*10^3, (z_grille)*10^3, IMAGE_env_norm_log);
% colormap(gray), caxis([-60 0]); 
% c = colorbar;
% c.Label.String = 'dB';
% c.Label.Position=[0.5 -60.8 0];
% c.Label.Rotation = 0;
% xlabel('Axe latéral de la sonde [mm]'); ylabel('Profondeur [mm]'); 
% title('Image mode B');
    