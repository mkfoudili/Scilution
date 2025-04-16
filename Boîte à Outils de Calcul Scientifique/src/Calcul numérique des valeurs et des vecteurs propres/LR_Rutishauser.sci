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

//Fonction de LR et Rutishauser
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

    clf;
    plot(1:N, convergence, '-o');
    xlabel("N");
    ylabel("Convergence de la methode du LR");
    title("Convergence de la methode du LR");
    xtitle("Convergence de la methide du LR");
endfunction