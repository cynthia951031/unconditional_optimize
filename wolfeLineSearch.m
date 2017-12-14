function alpha = wolfeLineSearch(objFunc,objFuncValue,x,dx,dir)
    alphaMax     = 1;
    alpha        = alphaMax;
    fac          = 1/2; 
    c_1          = 1e-1;
    c_2          = 2e-1;


    while objFunc(x+alpha*dir) > objFuncValue + c_1*alpha*dir'*dx && abs(numDiff(objFunc, x+alpha*dir)'*dir) < -c_2*dx'*dir
        alpha = fac*alpha;

        if alpha < 10*eps
            error('Error in Line search - alpha close to working precision');
        end

    end
end