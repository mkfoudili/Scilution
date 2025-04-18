//Fonctions de la sécante
function rNew = U(rOld , rOlder) // rn+1 = rn - f(rn)(rn - rn-1)/(f(rn) - f(rn-1))
        
    rNew = rOld - f(rOld)*(rOld - rOlder)/(f(rOld) - f(rOlder))
    
endfunction

function rn = SECANTE( f , r0 , r1 , n )
    
    rOlder = r0
    rOld = r1
    rNew = U(rOld , rOlder)
    
    mprintf("the new rn : %f\n" , rNew)
    
    for i=1:n-1
        
        rOlder = rOld
        rOld = rNew
        rNew = U(rOld , rOlder)
        
        mprintf("the new rn : %f\n" , rNew)    
        
        if isnan(rNew) then
            rNew = rOld
            break
        end
        
    end
    
    rn = rNew