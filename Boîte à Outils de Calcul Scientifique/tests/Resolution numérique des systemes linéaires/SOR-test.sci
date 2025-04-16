A = [2 1 0; 1 4 1; 0 2 4]; //A
b = [1; 0; 0]; //b
e = 0.001;
w = 0.5;

if(det(A) == 0) then
    mprintf("Matrice non inversible");
else
    [T, N] = SOR(A, b, w , e);
    mprintf("T execution: %f s\nNbr iterations: %d\n\n",T,N);
end