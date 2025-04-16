// Methodes devloppées : Dichotomie, Newton, SECANTE
//*************************************************************************************

funcprot(0);

function y = f(x)
    //entrer la formule de la fonction ici
    y = x^2 + 2*x - 1; // racine dans [0,1]
    //decommenter pour la deuxieme methode Newton et la methode de la sécante
    //y = exp(-x) - x
    
endfunction

function y = df(x) //After choosing your function , write its derivative here (PS : there should be something in scilab that can make this process automatic , but I couldn't find it)

    //y = 2*x + 2
    y = -exp(-x) - 1

endfunction

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

//Fonctions de methode de Newton
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
endfunction


//**********************************************************************************************
//*******************************   PROGRAMME PRINCIPAL  ***************************************

[rn , n] = DICHOTOMIE(f , 0 , 1 , 0.005);
mprintf("la racine acquis est : %f after %d iterations" , rn , n)

//change la fonction f avant d'executer
[rn , n] = NEWTON(f , df , 1 , 0.005)
mprintf("la racine acquis est : %f after %d iterations" , rn , n)

//change la fonction f avant d'executer
rn = SECANTE(f , 0 , 1 , 10)
mprintf("la racine acquis est : %f " , rn)
