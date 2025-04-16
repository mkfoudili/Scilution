//Fonction de puissances inverses

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