function [bin_call_AM, bin_put_AM] = binomial_american(S0,K,r,sigma,T)

% MODELLED AFTER "american.m" IN HIGHAM, D.J.: "NINE WAYS TO IMPLEMENT THE
% BINOMIAL METHOD FOR OPTION VALUATION IN MATLAB", SIAM REVIEW, 44(4),
% 2002, 661-677.

tic
M = 5000;
dt = T/M;
A = 0.5*(exp(-r*dt)+exp((r+sigma^2)*dt)); 
u = A+sqrt(A^2-1);
d = 1/u;
p = (exp(r*dt)-d)/(u-d);

dpowers = d.^((M:-1:0)');
upowers = u.^((0:M)');
scale1 = p*exp(-r*dt);
scale2 = (1-p)*exp(-r*dt);

bin_call_AM = max(S0*dpowers.*upowers-K,0);
for i = M:-1:1
    Si = S0*dpowers(M-i+2:M+1).*upowers(1:i);
    bin_call_AM = max(max(Si-K,0),scale1*bin_call_AM(2:i+1) + scale2*bin_call_AM(1:i));
end

bin_put_AM = max(K-S0*dpowers.*upowers,0);
for i = M:-1:1
    Si = S0*dpowers(M-i+2:M+1).*upowers(1:i);
    bin_put_AM = max(max(K-Si,0),scale1*bin_put_AM(2:i+1) + scale2*bin_put_AM(1:i));
end
end