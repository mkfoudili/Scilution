getd('../src/');
funcprot(0);

function y = f(x)
    y = exp(-x) - x
    
endfunction

function y = df(x)
    y = -exp(-x) - 1

endfunction

[rn , n] = NEWTON(f , df , 1 , 0.005)
mprintf("la racine acquis est : %f after %d iterations" , rn , n)