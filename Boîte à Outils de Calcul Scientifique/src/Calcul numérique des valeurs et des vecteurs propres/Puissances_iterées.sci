//Fonction de puissances iterées

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
