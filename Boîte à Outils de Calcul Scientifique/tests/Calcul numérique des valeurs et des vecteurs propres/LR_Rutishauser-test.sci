getd('../src/');
funcprot(0);

M=[4    1;1    4]
[AN, LN, c] = LR(M, 10)
mprintf("La matrice obtenue apres l''application du LR sur La matrice: \n")
disp(AN)
mprintf("Les vecteurs propres : \n")
disp(LN)