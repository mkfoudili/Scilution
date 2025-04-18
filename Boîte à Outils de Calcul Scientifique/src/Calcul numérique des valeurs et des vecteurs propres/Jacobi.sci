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


//Fonction de Jacobi
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

    clf;
    plot(1:N, convergence, '-o');
    xlabel("N");
    ylabel("Convergence de la methide du Jacobi");
    title("Convergence de la methide du Jacobi");
    xtitle("Convergence de la methide du Jacobi");
endfunction