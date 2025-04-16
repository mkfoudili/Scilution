getd('../src/');
funcprot(0);

function y = f(x)
    y = x^2 + 2*x - 1; // racine dans [0,1]
    
endfunction

function y = df(x)
    y = 2*x + 2
endfunction

[rn , n] = DICHOTOMIE(f , 0 , 1 , 0.005);
mprintf("la racine acquis est : %f after %d iterations" , rn , n)