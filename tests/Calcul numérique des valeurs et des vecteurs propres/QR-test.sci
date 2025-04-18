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