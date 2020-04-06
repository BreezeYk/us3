function x = BPSK1( OG , paramEm  , paramField  , nC , win , binary )

    paramEm.nS  = length( sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) ) ) 
    x = zeros( size(OG,1) , paramEm.nS*size(OG,2)  ) ; 
    w=(tukeywin( length( sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( 1,1 ) )  , 0.4  ) );
  % zerop = 2*size( paramField.rep_impuls , 2 ) + nC*paramEm.nS ;
    zerop = nC*paramEm.nS*paramEm.NbCycles ; 
   % zerop = 1600 ;
    c = floor( paramEm.nS*size(OG,2)  / zerop ) ;
    for j = 1 : size(OG,1)
      %  if j >= size( OG , 1 )
      %  else
      %      f= 1 ;
      %  for i = 1 : c
        for i = 1 : size(OG,2)
            
           % subC = zeros(nC*paramEm.nS,1) ; 
              %  for k = 1:nC
                %  if f >= size( OG , 1 )
                 %       break
                 % else
                 %   if k == 1
                 %       subC(k:k+paramEm.nS) = [ sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,f ) ] ;
                 %   else
                  %      subC( (k-1)*paramEm.nS : (k-1)*paramEm.nS + paramEm.nS) = [ sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,f ) ] ;
                  %  end
                 %  f = f + 1 ;
                %  end
              %  end
            %  if i >= size(OG,1)
            %      break
            %  else
            if binary == 1
                
                if i == 1 
                    x( j , i : i + paramEm.nS + zerop - 1  ) = [OG( j,i ) zeros(zerop+size( sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) ) , 2) -1  , 1)'] ;
                else
                    x( j , (i-1)*( paramEm.nS +zerop) : (i-1)*(paramEm.nS+zerop) + paramEm.nS + zerop - 1  ) = [OG( j,i ) zeros(zerop + size( sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) ) , 2) -1  ,1)']  ;
                end
                
            else
                    if i == 1 

                        if win == 1
                             x( j , i : i + paramEm.nS + zerop - 1 ) = [w'.*sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,i ) zeros(zerop,1)'] ;
                        else
                             x( j , i : i + paramEm.nS + zerop - 1  ) = [sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,i ) zeros(zerop,1)'] ;
                        end
                    else
                        if win ==1
                            x( j , (i-1)*( paramEm.nS +zerop) : (i-1)*(paramEm.nS+zerop) + paramEm.nS + zerop - 1  ) = [w'.*sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,i ) zeros(zerop,1)']  ;
                        else
                            x( j , (i-1)*( paramEm.nS +zerop) : (i-1)*(paramEm.nS+zerop) + paramEm.nS + zerop - 1  ) = [sin( 2*pi*paramField.f0*(0:1/paramField.fs: paramEm.NbCycles/paramField.f0) )*OG( j,i ) zeros(zerop,1)']  ;
                        end
                    end

                end
        end
    end
end