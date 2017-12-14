function [x, fval] = augmented_lagrangian_method(f, ce, ci, x0)
    prec          = 10e-8;
    x             = x0;
    numOfIter     = 100;
    lambda        = ones(length(ce(x)) + length(ci(x)), 1); %lambda(i) CORRESPOND TO CONDITIONS
    penalty       = 0.1;
    penaltyFactor = 10;
    method        = 'fminunc';

    eta           = lambda / penalty;
    neq           = length(ce); % OBTAIN THE QUANTITY OF EUQAL AND INEQUAL CONDITIONS
    %LAGRANGIAN FUNCTION
    if ~isempty(ci(x))
        L = @(x, lambda, penalty, eta) (f(x) - ce(x) * lambda(1 : neq) +...
            0.5 * penalty * norm(max(ce(x), 0))^2 + ...
            sum(0.5 * penalty * (min(ci(x) - eta(neq + 1: end)', 0).^ 2) ...
            - eta(neq + 1: end)' .^ 2));
    else
        L = @(x, lambda, penalty, eta) (f(x) - ce(x) * lambda(1 : neq) +...
            0.5 * penalty * norm(max(ce(x), 0))^2);
    end

    for iter = 1: numOfIter
        %OBJECT FUNCTION
        objFun = @(x) L(x, lambda, penalty, eta);
        %UNCONDITIONAL OPTIMIZE
        xold = x;
        [x, grad] = optim(objFun, xold, method);
        %UPDATE LAMBDA
        lambda(1: neq) = (-penalty * min(ce(x) - (eta(1 : neq))', 0))';
        if ~isempty(ci(x))
            lambda(neq + 1: end) = (-penalty * min(ci(x) - (eta(neq + 1: end))', 0))';
        end
        %UPDATE PENALTY
        penalty = penalty * penaltyFactor;
        %NORM VERIFICATION
        normNum = norm(ce(x))^2;
        if ~isempty(ci(x))
            normNum = normNum + norm(min(ci(x), eta(neq+1:end)'))^2;
        end
        
        if sqrt(normNum) <= prec
            break;
        end
        
        %{
        %KKT VERIFICATION
        DxL = @(s, lambda, rho)(dF(s) - dC(s)*lambda + rho*dC(s)*C(s)');
        DxPr = @(s, lambda, rho)(almProj(s, DxL(s, lambda, rho), lb, ub));
        P = @(s, lambda, rho)deal(L(s, lambda, rho), s-DxPr(s, lambda, rho));
        if neq && nineq
            KKT = @(s, lambda, rho)deal(s - DxPr(s, lambda, rho), ce(s(1:nvar)), ci(s(1:nvar)));
        elseif nineq
            KKT = @(s, lambda, rho)deal(s - DxPr(s, lambda, rho), 0.0, ci(s(1:nvar)));
        elseif neq
            KKT = @(s, lambda, rho)deal(DxL(s, lambda, rho), ce(s(1:nvar)), 0.0);
        else
            KKT = @(s, lambda, rho)deal(df(s(1:nvar)), 0.0, 0.0);
        end
        %}
        %ETA UPDATE
        eta = lambda / penalty;
        fprintf(1,'\nALM Iteration %d: R=%f, F=%f\n',iter, penalty, f(x));
    end
    fval = f(x);
end