function [xnew, grad] = optim(f, xold, method)
    switch lower(method)
        case 'bfgs'
            xnew = bfgs(xold, f);
            grad = numDiff(f,xnew);
        case 'fminunc'
            [xnew, ~, ~, ~, grad] = fminunc(f, xold);
    end
end