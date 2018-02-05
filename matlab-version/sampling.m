function [z,q,G] = sampling(m,t,gamma,type,strategy)

% THIS FUNCTION INCLUDES THE SAMPLING STRATEGIES FOR {z(i),q(i)}

if strcmp(type,'Curran')
    
    x = 1:m; x=x';
    z = norminv((x-0.5)/m);
    q = 1/m*ones(m,1);

elseif strcmp(type,'new sampling strategy')
    
    x = 1:m/2; x=x';
    q = (x-0.5).^gamma/m;

    if mod(m,2) == 0        
        q = [q; flip(q)];        
    else        
        q_mid = (ceil(m/2)-0.5).^gamma/m;
        q = [q; q_mid; flip(q)];        
    end
    
    q = q./sum(q);

    Z = [-Inf; norminv(cumsum(q(1:end-1))); Inf];
    
    options = optimoptions('fmincon','Display','final','TolX',1e-15,...
        'TolFun',1e-15,'TolCon',1e-8,'MaxIter',1e4,'MaxFunEvals',1e5);

    if strcmp(strategy,'OPT strategy')
        
        tic    
        z = fmincon(@(z) (q'*z.^4-3).^2,zeros(m,1),[],[],[],[],...
            Z(1:end-1),Z(2:end),@(z) nonlcon(z,q),options);        
        toc
    
    elseif strcmp(strategy,'ADJ strategy')
        tic
        if mod(m,2) == 0
           z = norminv([q(1)/2; cumsum(q(1:end-1))+q(2:end)/2]);
           a = (3/2 - q(1:m/2)'*z(1:m/2).^4)/...
               ((z(1)^2-z(2)^2)*(z(1)^2-z(m/2)^2));
           b = (3/2 - q(1:m/2)'*z(1:m/2).^4)/...
               ((z(2)^2-z(1)^2)*(z(2)^2-z(m/2)^2));
           q_hat = [q(1)+a; q(2)+b; q(3:m/2-1); q(m/2)-a-b];
           q = [q_hat; flip(q_hat)];
        else
        end
        toc
    
    elseif strcmp(strategy,'FPM strategy')        
        tic        
        z = fmincon(@(z) sum(abs((q'*bsxfun(@(z,Z) max(z-Z,0),z,...
            Z(2:end-1)'))-(1/(sqrt(2*pi))*exp(-Z(2:end-1).^2/2)-...
            Z(2:end-1).*(1-normcdf(Z(2:end-1))))')),zeros(m,1),...
            [],[],[],[],Z(1:end-1),Z(2:end),@(z) nonlcon(z,q),options);        
        toc        
    end        
end

G = bsxfun(@(t,z) t.*z, sqrt(t),z);

end