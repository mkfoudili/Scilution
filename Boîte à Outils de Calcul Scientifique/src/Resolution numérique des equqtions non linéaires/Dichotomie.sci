//Fonction de la dichotomie
function [rn , n] = DICHOTOMIE( f , a , b , e )
    cptMAX = 500
    n = 0
    
    while abs(b - a)/2 > e && n < cptMAX do        
        rn = (b + a)/2        
        mprintf("[%f , %f] and rn = %f \n", a , b , rn)        
        fa = f(a)
        fb = f(b)
        fr = f(rn)       
        if(fa*fr < 0) then b = rn ;
            else a = rn ;
        end 
        n = n + 1
        
    end
    
    if abs(b - a)/2 > e || f(a)*f(b) > 0 then mprintf("le programme est arreté aprés %d iterations \n\n" , n) end
    
endfunction