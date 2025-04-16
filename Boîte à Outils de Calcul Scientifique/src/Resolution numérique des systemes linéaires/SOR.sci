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

//Fonction de SOR
function [T,N]= SOR(A,b,w,e) 
    
    if(w <= 0 || w >= 2) then
        mprintf("change w!");
        T = 0;
        N = 0;
        exit();
    end
    
    tic();
    
    cpt = 100;
    N = 0;
    
    Lsor = inv( inv(w)*D(A) - E(A) ) * ( F(A) + (1 - inv(w))*D(A) ); //LSOR: inv(M)*N = inv(1/w * D - E)*(F + (1 - 1/w)D)
    C = inv( inv(w)*D(A) - E(A) ) * b; //inv(M)*b = inv(1/w * D - E)*b
    
    Xold = [1; 1; 1]; //X0
    Xnew = Lsor * Xold + C; //X1
    
    erreur1 = norm(Xnew - Xold, 1) * norm(Lsor , 1) / (1 - norm(Lsor , 1)) //erreur 1
    erreurinfini = norm(Xnew - Xold, "inf") * norm(Lsor , "inf") / (1 - norm(Lsor , "inf")) //erreur infini
    
    //Boucle principale
    while erreur1 > e && erreurinfini > e && N < cpt 
        
        Xold = Xnew; //Xn+1
        Xnew = Lsor * Xold + C; //calcul Xn+1
        
        erreur1 = norm(Xnew - Xold, 1) * norm(Lsor , 1) / (1 - norm(Lsor , 1))
        erreurinfini = norm(Xnew - Xold, "inf") * norm(Lsor , "inf") / (1 - norm(Lsor , "inf")) 
        N = N + 1;
        
    end
    
    T = toc();//Fin
    
    if(N == cpt) then mprintf("Methode peut etre divergente");
    else
        mprintf("la valeur X approché trouvee par SOR avec w=%f est :\n",w);
        disp(Xnew);
    end
    
endfunction