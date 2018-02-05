function [bin_call, bin_put] = binomial(S0,K,r,sigma,T,t,steps)

% MODELLED AFTER "EURO5.m" IN HIGHAM, D.J.: "NINE WAYS TO IMPLEMENT THE
% BINOMIAL METHOD FOR OPTION VALUATION IN MATLAB", SIAM REVIEW, 44(4),
% 2002, 661-677.

tic

if strcmp(steps,'m')
M = (length(t)-2);
elseif strcmp(steps,'5000')
M = 5000;
end
    
dt = T/M;
A = 0.5*(exp(-r*dt)+exp((r+sigma^2)*dt)); 
u = A+sqrt(A^2-1);
d = 1/u;
p = (exp(r*dt)-d)/(u-d);
q = 1-p;
bin_call = max(S0*d.^((M:-1:0)').*u.^((0:M)')-K,0);
bin_put = max(K-S0*d.^((M:-1:0)').*u.^((0:M)'),0);
for i = M:-1:1
    bin_call = p*bin_call(2:i+1) + q*bin_call(1:i);
end
bin_call = exp(-r*T)*bin_call;
for i = M:-1:1
    bin_put = p*bin_put(2:i+1) + q*bin_put(1:i);
end
bin_put = exp(-r*T)*bin_put;
toc
end