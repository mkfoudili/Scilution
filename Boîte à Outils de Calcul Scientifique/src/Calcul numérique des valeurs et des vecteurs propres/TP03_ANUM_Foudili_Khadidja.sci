// Methodes devloppées : Méthodes partielles (la puissance iterées, inverse), Méthodes itératives globales (Jacobi, LR et Rutishauser, QR)

//*************************************************************************************

funcprot(0); 
//Fonctions de generation de A
function A = genererA(d)
    A = zeros(d, d);
    for i = 1:d
        for j = 1:d
            A(i, j) = 1 / (i + j - 1);
        end
    end 
endfunction

function B = genererB(d)
//Fonctions de generation de B
    B = zeros(d, d);
    for i = 1:d
        B(i, i) = 2;
        if( i < d ) then B(i, i+1) = -1; end
        if( i > 1 ) then B(i, i-1) = -1; end
    end
endfunction

// Methode de la puissance iterée

function [vp, vect, error] = PuissanceIteree(A, x0, N)
    xOld = x0 / norm(x0);
    vpOld = 0;
    vp = 0;
    vp_evolution = [];

    for k = 1:N
        xNew = A * xOld;
        vp = (xOld' * xNew) ;
        error = abs(vp - vpOld);
        vp_evolution = [vp_evolution, vp];
        xNew = xNew / norm(xNew);       
        vpOld = vp;
        xOld = xNew
    end
    vect = xNew;
    
    plot(1:N, vp_evolution, '-o');
    xlabel("N");
    ylabel("Plus Grande Valeur Propre");
    title("VP evolution par iteration");
    xtitle("Methode de la puissance iterée");
endfunction

// Methode de la puissance inverse
function [vp, vect, error] = PuissanceInverse(A, x0, N)
    xOld = x0 / norm(x0);
    vpOld = 0;
    vp = 0;
    vp_evolution = [];
    A = inv(A)

    for k = 1:N
        xNew = A * xOld;
        vp = (xOld' * xNew) ;
        error = abs(vp - vpOld);
        vp_evolution = [vp_evolution, vp];
        xNew = xNew / norm(xNew);       
        vpOld = vp;
        xOld = xNew
    end
    vect = xNew;
    vp = 1/vp
    //decommenter pour tracer
    /*
    plot(1:N, vp_evolution, '-o');
    xlabel("N");
    ylabel("Plus Grande Valeur Propre");
    title("VP evolution par iteration");
    xtitle("Methode de la puissance iterée"); */
endfunction

// Methode de Jacobi
function [AN, QN, convergence] = Jacobi(A, N)    
    d = size(A, 1);
    QN = eye(d, d);//multiplication des matrices (ses colonnes sont les vecteurs propres associés)
    AN = A;
    convergence = [];
    for k = 1:N
        [max_val, p, q] = MaxHorsDiagonal(AN);
        convergence = [convergence, max_val];
        Q = construireQPQ(AN, p, q)
        AN = Q' * AN * Q; //la transpose est l'inverse
        QN = QN * Q;
    end
     
    //decommenter pour tracer  
    /*
    clf;
    plot(1:N, convergence, '-o');
    xlabel("N");
    ylabel("Convergence de la methide du Jacobi");
    title("Convergence de la methide du Jacobi");
    xtitle("Convergence de la methide du Jacobi");
    */
endfunction

function [max_val, p, q] = MaxHorsDiagonal(A)
    d = size(A, 1);
    max_val = 0;
    p = 0; q = 0;
    for i = 1:d-1 
        for j = i+1:d
            if abs(A(i, j)) > max_val
                max_val = abs(A(i, j));
                p = i; q = j;
            end
        end
    end
endfunction

function J = construireQPQ(A, p, q)
    d = size(A, 1);
    if A(p, p) == A(q, q)
        cosTeta = 1 / sqrt(2);
        sinTeta = 1 / sqrt(2);
    else
        b = (A(q, q) - A(p, p)) / (q * A(p, q));
        t = sign(b) / (abs(b) + sqrt(1 + b^2));
        cosTeta = 1 / sqrt(1 + t^2);
        sinTeta = t * cosTeta;
    end
    J = eye(d, d);
    J(p, p) = cosTeta; J(q, q) = cosTeta;
    J(p, q) = sinTeta; J(q, p) = -sinTeta;
endfunction

// Methode de QR
function [AN, QN] = QR(A, N)
    d = size(A, 1);
    QN = eye(d, d);
    AN = A;
    for k = 1:N
        [Q, R] = qr(AN);
        //la fonction de decomposition existe deja sur scilab
        AN = R * Q;
        QN = QN * Q;
    end
endfunction

// Methode de LR et Rutishauser
function [AN, LN, convergence] = LR(A, N)
    d = size(A, 1);
    LN = eye(d, d);
    AN = A;
    convergence = [];

    for k = 1:N
        [L, U] = lu(AN);
         //la fonction de decomposition existe deja sur scilab
        AN = U * L;
        LN = LN * L;
        triangInf = tril(AN, -1);   // Extract strictly lower triangular part
        maxHorsDiagonal = max(abs(triangInf));
        convergence = [convergence, maxHorsDiagonal];
    end
    
    //decommenter pour tracer  
    /*
    clf;
    plot(1:N, convergence, '-o');
    xlabel("N");
    ylabel("Convergence de la methode du LR");
    title("Convergence de la methode du LR");
    xtitle("Convergence de la methide du LR");
    */
endfunction

//**********************************************************************************************
//*******************************   PROGRAMME PRINCIPAL  ***************************************
//soit A une mtrice:
A = genererA(10);
B = genererB(10);

mprintf("\n\nMatrice A :\n");
disp(A);

mprintf("\n\nMatrice A :\n");
disp(B);

x0 = diag(A) 
[vp, X, e] = PuissanceIteree(A, x0, 12)

mprintf("La plus grande valeur propre de A est : %f \n", vp)
mprintf("erreur: %f\n",e)
mprintf("Son vecteur propre X :  \n");
disp(X)

x0 = diag(B) 
[vp, X, e] = PuissanceIteree(B, x0, 12)

mprintf("La plus grande valeur propre de B est : %f \n", vp)
mprintf("erreur: %f\n",e)
mprintf("Son vecteur propre X :  \n");
disp(X)

x0 = diag(A) 
[vp, X, e] = PuissanceInverse(A, x0, 12)

mprintf("La plus petite valeur propre de A est : %f \n", vp)
mprintf("erreur: %f\n",e)
mprintf("Son vecteur propre X :  \n");
disp(X)

x0 = diag(B) 
[vp, X, e] = PuissanceInverse(B, x0, 12)

mprintf("La plus petite valeur propre de B est : %f \n", vp)
mprintf("erreur: %f\n",e)
mprintf("Son vecteur propre X :  \n");
disp(X)

[AN, QN, c] = Jacobi(A, 12)
mprintf("La matrice obtenue apres l''application de jacobi sur A : \n")
disp(AN)
mprintf("Les vecteurs propre : \n")
disp(QN)

[AN, QN, c] = Jacobi(B, 12)
mprintf("La matrice obtenue apres l''application de jacobi sur A : \n")
disp(AN)
mprintf("Les vecteurs propres : \n")
disp(QN)

[AN, QN] = QR(A, 10)
mprintf("La matrice obtenue apres l''application du QR sur A : \n")
disp(AN)
mprintf("Les vecteurs propres : \n")
disp(QN)

[AN, QN] = QR(B, 10)
mprintf("La matrice obtenue apres l''application du QR sur B : \n")
disp(AN)
mprintf("Les vecteurs propres : \n")
disp(QN)

M=[4    1;1    4]
[AN, LN, c] = LR(M, 10)
mprintf("La matrice obtenue apres l''application du LR sur La matrice de l''exercice 4 : \n")
disp(AN)
mprintf("Les vecteurs propres : \n")
disp(LN)
