function y = BSPK( OG , paramEm , paramField , L , win )

    exArray = zeros( size(OG,1) , paramEm.nS*L) ;
    
    
    w = tukeywin( length( paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0)) ) , 0.8 ) ;
    
    for i = 1:size(OG,1)
       for j = 1 : L
           
           if j == 1
               if  win == 1
                   exArray( i , j:j+paramEm.nS-1 ) = w'.*paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0))*OG(i,j) ;
               else
                   exArray( i , j:j+paramEm.nS-2 ) = paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0))*OG(i,j) ;
               end
           else
                   if win == 1
                    exArray( i , (j-1)*paramEm.nS : (j-1)*paramEm.nS+paramEm.nS-1 ) = w'.*paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0))*OG(i,j) ;
                   else
                       exArray( i , (j-1)*paramEm.nS : (j-1)*paramEm.nS+paramEm.nS-2 ) = paramEm.excitation(1:floor(paramEm.NbCycles*paramField.fs/paramField.f0))*OG(i,j) ;
                   end
                   end
        end
       
    end
    
    y = exArray ;
    
end