funcprot(0); 
function result = D(A) 
    
    d = diag(A);
    result = diag(d);
    
endfunction

function result = E(A) 
    
    L = tril(A); 
    result = -(L - D(A)); 
    
endfunction

function result = F(A) 
    
    U = triu(A); 
    result = -(U - D(A)); 
    
endfunction

//Fonction de Gauss Seidel
function [T,N]= GaussSeidel(A,b,e) 
    
    tic(); 
    
    cpt = 100; 
    N = 0;
    
    Lgs = inv( D(A) - E(A) ) * F(A); //LGS: inv(M)*N = inv(D - E)*F
    C = inv( D(A) - E(A) ) * b; //inv(M)*b = inv(D - E)*b
    
    Xold = [1; 1; 1]; //X0
    Xnew = Lgs * Xold + C; //X1
    
    erreur1 = norm(Xnew - Xold, 1) * norm(Lgs , 1) / (1 - norm(Lgs , 1)) //erreur 1
    erreurinfini = norm(Xnew - Xold, "inf") * norm(Lgs , "inf") / (1 - norm(Lgs , "inf")) //erreur infini
    
    //Boucle principale
    while erreur1 > e && erreurinfini > e && N < cpt 
        Xold = Xnew; //Xn+1
        Xnew = Lgs * Xold + C; //Xn+1
        
        erreur1 = norm(Xnew - Xold, 1) * norm(Lgs , 1) / (1 - norm(Lgs , 1)) 
        erreurinfini = norm(Xnew - Xold, "inf") * norm(Lgs , "inf") / (1 - norm(Lgs , "inf")) 
        N = N + 1;
        
    end
    
    T = toc(); //Fin
    
    if(N == cpt) then mprintf("Methode peut etre divergente");
    else
        mprintf("X:\n");
        disp(Xnew);
    end
    
endfunction