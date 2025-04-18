getd('../src/');
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

//soit A une mtrice:
A = genererA(10);
B = genererB(10);

mprintf("\n\nMatrice A :\n");
disp(A);

mprintf("\n\nMatrice A :\n");
disp(B);

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