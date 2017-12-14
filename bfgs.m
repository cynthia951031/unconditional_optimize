function x = bfgs(init_x, f)
x       = init_x'; 
H       = eye(length(init_x));  

objFunc         = f;
objFuncValue    = objFunc(x);
oldObjFuncValue = objFuncValue + 1;
dx = numDiff(objFunc,x);

iter      = 0;
numOfIter = 100;
prec      = 1e-10;

while iter < numOfIter && abs((oldObjFuncValue-objFuncValue)/objFuncValue)>prec && norm(dx)>prec
    iter = iter + 1;
    oldObjFuncValue = objFuncValue;
    dir = -(H\dx); 
    alpha = wolfeLineSearch(objFunc,objFuncValue,x,dx,dir);
    
    p = alpha*dir;
    x = x + p;
    objFuncValue = objFunc(x);
    dx_old = dx;
    dx = numDiff(objFunc,x);
    q = dx-dx_old;

    H = H + (1+(p'*H*p)/(p'*q))*(q*q')/(p'*q) - (q*p'*H+H*p*q')/(p'*q);
    
    %fprintf(1,'Iteration %d: alpha_min=%f, OF=%f\n', iter, alpha, objFuncValue);
    
end
%fprintf(['\n' num2str(iter) ' iteration(s) performed to converge\n'])
%fprintf(1,'Final solution: \n');
%display(x);
%display(objFunc(x));
end