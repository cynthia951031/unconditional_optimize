function dp = numDiff(func,p)
    h  = 1e-8;
    dp = NaN*p;

    oldObjFuncValue = func(p);

    for i = 1:numel(p)

        p_new    = p;
        p_new(i) = p_new(i) + h;

        newObjFuncValue = func(p_new);

        dp(i) = (newObjFuncValue-oldObjFuncValue)/h;

    end
end