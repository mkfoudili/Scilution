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

//Fonction de Jaconbi
function [T,N]= Jacobi(A,b,e) 
    
    tic(); 
    
    cpt = 100;
    N = 0;
    
    Lj = inv(D(A)) * (E(A) + F(A)); //LJ: inv(M)*N = inv(D)*(E+F)
    C = inv(D(A)) * b; //inv(M)*b = inv(D)*b
    
    Xold = [1; 1; 1]; //X0
    Xnew = Lj * Xold + C; //X1
    
    erreur1 = norm(Xnew - Xold, 1) * norm(Lj , 1) / (1 - norm(Lj , 1)) //erreur 1
    erreurinfini = norm(Xnew - Xold, "inf") * norm(Lj , "inf") / (1 - norm(Lj , "inf")) //erreur infini
    
    //Boucle principale
    while erreur1 > e && erreurinfini > e && N < cpt 
        
        Xold = Xnew; //Xn+1
        Xnew = Lj * Xold + C; //Xn+1
        
        erreur1 = norm(Xnew - Xold, 1) * norm(Lj , 1) / (1 - norm(Lj , 1)) //maj erreur 1
        erreurinfini = norm(Xnew - Xold, "inf") * norm(Lj , "inf") / (1 - norm(Lj , "inf")) //maj erreur infini
        N = N + 1;
        
    end
    
    T = toc(); //Fin
    
    if(N == cpt) then mprintf("methode peut etre divergente");
    else
        mprintf("la valeur X approché trouvee par Jacobi est :\n");
        disp(Xnew);
    end
    
endfunction