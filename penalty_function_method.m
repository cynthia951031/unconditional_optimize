function [x, fval] = penalty_function_method(f, ce, ci, x0)
    prec           = 1;
    numOfIter      = 100;
    penaltyFactor  = 10;
    penalty        = 0.1;
    x              = x0;
    method         = 'fminunc';
    
    for iter = 1: numOfIter
        %OBJECT FUNCTION
        objFun = @(x) (f(x) + 0.5 * penalty * (norm(max(ce(x), 0))^2 + norm(min(ci(x), 0))^2));
        %UNCONDITIONAL OPTIMIZE
        xold = x;
        [x, grad] = optim(objFun, xold, method);
        %UPDATE PENALTY
        penalty = penalty * penaltyFactor;
        %NORM VERIFICATION
        if norm([ce(x), ci(x)]) < prec
            break;
        end

        fprintf(1,'\nPenalty Iteration %d: R=%f, F=%f, Err=%f\n',iter, penalty, f(x), norm([ce(x), ci(x)]));
    end
    fval = f(x);
end



