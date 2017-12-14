function [x, fval] = barrier_function_method(f, ci, x0)
    prec           = 1e-7;
    numOfIter      = 100;
    penaltyFactor  = 0.5;
    penalty        = 1;
    x              = x0;
    method         = 'fminunc';
    barrier_func   = 'inverse'; %'inverse' / 'log'
    
    for iter = 1: numOfIter
        switch lower(barrier_func)
            case 'inverse'
                objFun = @(x) f(x) + penalty * sum(ci(x) .^ -1);
            case 'log'
                objFun = @(x) f(x) - penalty * sum(log(ci(x)));
        end
        %UNCONDITIONAL OPTMIZE
        xold = x;
        [x, grad] = optim(objFun, xold, method);
        if penalty * sum(ci(x) .^ -1) <= prec
            break;
        end
        %UPDATE PENALTY
        penalty = penalty * penaltyFactor;
        fprintf(1,'\nPenalty Iteration %d: R=%f, F=%f, Err=%f\n',iter, penalty, f(x), penalty * sum(ci(x) .^ -1));
    end
    fval = f(x);
end