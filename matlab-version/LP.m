function [P,fval,exitflag,h] = LP(m,z,q,t,gamma,type)

% THIS FUNCTION INCLUDES THE LINEAR PROGRAMMING (LP) PROBLEM.

if strcmp(type,'new sampling strategy')
    
    u = ones(m,1);
    r = z.^2;
    h = t(3:end) - t(2:end-1);
    alpha = h(1:end)./t(2:end-1);
    beta = 1./sqrt(1+alpha);

    temp = bsxfun(@(m,z) m.*z, ones(m),z);
    F  = cell(1,numel(h)); P = F;

    for i=1:numel(h)
        
        F{i} = (abs(temp' - beta(i)*temp)).^3;
        F{i} = F{i}'; F{i} = F{i}(:);
        
    end

    Aeq = [kron(eye(m),u');
           kron(eye(m),z');
           kron(eye(m),r');
           kron(q',eye(m))];
    
    Aeq = sparse(Aeq);

    options = optimoptions('linprog','Algorithm','dual-simplex',...
        'Display','final','TolCon',1e-9);
    
    tic
   
    for i = 1:numel(h)
        
        beq = [u;
               beta(i)*z;
               beta(i)^2*r+(1-beta(i)^2)*u;
               q];
           
           if gamma == 0
               
               [P{i},fval,exitflag] = linprog(F{i},[],[],Aeq,beq,...
                   zeros(m*m,1),[],[],options);
           
           else
               
               [P{i},fval,exitflag] = linprog(kron(q,ones(m,1)).*F{i},...
                   [],[],Aeq,beq,zeros(m*m,1),[],[],options);
           
           end
    
    P{i} = reshape(P{i},m,m); P{i} = P{i}';
    
    end
    
    toc
    
end

end