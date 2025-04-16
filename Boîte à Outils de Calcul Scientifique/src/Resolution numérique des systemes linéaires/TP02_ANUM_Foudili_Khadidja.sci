//Methodes devloppées : Jacobi, Gauss Seidel, et SOR

//*************************************************************************************
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

//**************************************************************************************

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

//**********************************************************************************************
//*******************************   PROGRAMME PRINCIPAL  ***************************************

A = [2 1 0; 1 4 1; 0 2 4]; //A
b = [1; 0; 0]; //b
e = 0.001;
w = 0.5;//pour la methode SOR

if(det(A) == 0) then
    mprintf("Matrice non inversible");
else
    //La valeur approchée pour les trois methodes     
    [T, N] = Jacobi(A, b, e);
    mprintf("T execution: %f s\nNbr iterations: %d\n\n",T,N);
    [T, N] = GaussSeidel(A, b, e);
    mprintf("T execution: %f s\nNbr iterations: %d\n\n",T,N);
    [T, N] = SOR(A, b, w , e);
    mprintf("T execution: %f s\nNbr iterations: %d\n\n",T,N);
end

