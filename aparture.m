function [emit_aperture , receive_aperture] = aparture( paramField )        

        emit_aperture = xdc_linear_array (paramField.N_elements_emission,...
                                          paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
                                          1, 1,paramField.focus);    % le focus qui est mis ici est temporaire et ne va pas servir a la focalisation, 
                                                                     % car apres en imposant les delays (signaux d'excitation) on supprime l'effet de ce focus


        receive_aperture = xdc_linear_array (paramField.N_elements_reception,...
                                          paramField.largeur_element, paramField.hauteur_element, paramField.kerf,...
                                          1, 1,paramField.focus);    % le focus qui est mis ici est temporaire et ne va pas servir a la focalisation, 
                                                                     % car apres en imposant les delays (signaux d'excitation) on supprime l'effet de ce focus


        % % 

        paramField.probe_specs = xdc_get(emit_aperture , 'rect'); % on recupere tous les parameters geometriques de la probe
        paramField.position_elements=paramField.probe_specs(8,:); % on recupere les coordonees des centres des elements de la probe
        
        xdc_impulse(emit_aperture, paramField.rep_impuls);      % en emission
        xdc_impulse(receive_aperture, paramField.rep_impuls);   % en reception

        % On supprime le focus en emission
        xdc_center_focus(emit_aperture,[0,0,0]);
        xdc_center_focus(receive_aperture,[0,0,0]);

        % On force les delays d'emission a 0 dans un premier temps
        xdc_focus_times( emit_aperture , 0, zeros(1,paramField.N_elements_emission));        % en emission 
        xdc_focus_times( emit_aperture , 0, zeros(1,paramField.N_elements_reception));    % en reception

        % On supprime le focus en emission
        xdc_center_focus(emit_aperture,[0,0,0]);

        xdc_apodization( emit_aperture,0, ones(1,paramField.N_elements_emission) );
        xdc_apodization( receive_aperture,0, ones(1,paramField.N_elements_reception) );