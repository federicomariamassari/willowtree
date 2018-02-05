function [S, willow_call, willow_put] = pricing(S0,K,r,sigma,m,q,G,P,t,T,strategy)

% PRICING FUNCTION FOR THE SAMPLING STRATEGIES

tic
if strcmp(strategy,'iCUR') || strcmp(strategy,'OPT') || strcmp(strategy,'FPM')
    temp = P{1};
    for i = 1:length(P)
        temp = temp*P{i};
        cond_exp = q'*temp;
    end
elseif strcmp(strategy,'CUR') || strcmp(strategy,'ADJ')
    cond_exp = q';
end

time = kron(t,ones(m,1));
S = S0*exp((r-0.5*sigma^2)*time+sigma*G);

X_call = max(S-K,0);
willow_call = exp(-r*T)*cond_exp*X_call(:,end);
   
X_put = max(K-S,0);
willow_put = exp(-r*T)*cond_exp*X_put(:,end);
toc       
end

