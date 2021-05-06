% ------------------------------------------------------------------------
% FINDING ZERO OF FUNCTION WITHIN A GIVEN BOUND BY BISECTION METHOD
% ------------------------------------------------------------------------

function [xzero,err] = m_Bisection(func,lb,ub)

errTol = 1e-6;

% Check if a zero exists in bound, i.e. func changes sign from lb to ub:
if abs(func(lb)) < errTol
    xzero = lb;
    err = abs(func(lb));
elseif abs(func(ub)) < errTol
    xzero = ub;
    err = abs(func(ub));
elseif func(lb)*func(ub) > 0
    error('Probably no zero within bound specified!')
else
    xzero = (lb + ub)/2;
    err = func(xzero);
    iter = 0;
    while iter < 1000
        iter = iter + 1;
        if func(lb)*func(xzero) < 0
            ub = xzero;
        else
            lb = xzero;
        end
        xzero = (lb + ub)/2;
        err = func(xzero);
        if abs(err) < errTol
            break
        end
    end
end

