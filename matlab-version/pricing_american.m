function [American_call, American_put] = pricing_american(m,q,t,G,P,S0,K,r,sigma,strategy)
tic
if strcmp(strategy,'iCUR') || strcmp(strategy,'OPT') || strcmp(strategy,'FPM')
    
    time = kron(t,ones(m,1));
    
    % AMERICAN CALL OPTION
    X_call = max(S0*exp((r-0.5*sigma^2)*time+sigma*G)-K,0);
    CV_call = NaN(m,numel(t)-2);
    CV_call(:,end) = max(X_call(:,end-1), exp(-r*(t(end)-t(end-1)))*P{end}*X_call(:,end));
    for i = size(CV_call,2)-1:-1:1
        CV_call(:,i) = max(X_call(:,i+1), exp(-r*(t(i+1)-t(i)))*P{i}*CV_call(:,i+1));
    end
    American_call = exp(-r*(t(2)-t(1)))*q'*CV_call(:,1);
    
    % AMERICAN PUT OPTION
    X_put = max(K-S0*exp((r-0.5*sigma^2)*time+sigma*G),0);
    CV_put = CV_call;
    CV_put(:,end) = max(X_put(:,end-1), exp(-r*(t(end)-t(end-1)))*P{end}*X_put(:,end));
    for i = size(CV_put,2)-1:-1:1
        CV_put(:,i) = max(X_put(:,i+1), exp(-r*(t(i+1)-t(i)))*P{i}*CV_put(:,i+1));
    end
    American_put = exp(-r*(t(2)-t(1)))*q'*CV_put(:,1);
toc
else
end