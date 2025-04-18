getd('../src/');
funcprot(0);

function y = f(x)
    y = exp(-x) - x
    
endfunction

function y = df(x)
    y = -exp(-x) - 1
endfunction

rn = SECANTE(f , 0 , 1 , 10)
mprintf("la racine acquis est : %f " , rn)