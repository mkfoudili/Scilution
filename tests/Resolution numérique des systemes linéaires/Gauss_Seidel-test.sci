getd('../src/');
funcprot(0);

A = [2 1 0; 1 4 1; 0 2 4]; //A
b = [1; 0; 0]; //b
e = 0.001;

if(det(A) == 0) then
    mprintf("Matrice non inversible");
else
    [T, N] = GaussSeidel(A, b, e);
    mprintf("T execution: %f s\nNbr iterations: %d\n\n",T,N);
end

