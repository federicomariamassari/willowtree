function [c,ceq] = nonlcon(z,q)
    
    % NONLINEAR EQUALITY CONTRAINTS - FMINCON INPUT
    
    c = [];
    ceq = [q'*z; q'*z.^2-1];
    
end