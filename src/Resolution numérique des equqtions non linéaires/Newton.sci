//Fonction de Newton
function rNew = g(rOld) // rn+1 = g(rn)
    
    rNew = rOld - f(rOld) / df(rOld) // rn+1 = rn + f(rn) / f'(rn)
    
endfunction

function [rn , n] = NEWTON( f , df , r0 , e )
    cptMAX = 50
    n = 1
    
    
    rOld = r0
    rNew = g(rOld)
    mprintf("the new rn : %f\n" , rNew)
    
    while abs(rNew - rOld) > e && n < cptMAX  do
        rOld = rNew
        rNew = g(rOld)        
        n = n + 1 ;   
        
        mprintf("the new rn : %f\n" , rNew)     
        
    end
    
    
    if abs(rNew - rOld) > e  then mprintf("le programme est arrete aprés %d iterations \n\n" , n) end
    
    rn = rNew
endfunction