
function ex = direcGene( fs , N )
    % y =  dirac(-0.000001:1/fs:0.000001) ;
    y =  dirac(0:1:N-1) ;
   idd = find( y == inf ) ;
   y( idd ) = 1*1e3 ;
    ex = y ;
end