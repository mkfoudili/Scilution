//Fonction de QR

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
