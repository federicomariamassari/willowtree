function [lookback_put, lookback_call] = lookback_american(m,t,S,q,P,r)

F0 = S(:,1);
F = cell(1,length(t)-1);
F_hat = F; [F_hat{:}] = deal(NaN(m,m));

F{1} = max(F0,S(:,2)); F_hat{1} = F{1};
F{1} = F{1}';
F{2} = bsxfun(@max,F{1},S(:,3)); F_hat{2} = F{2};

for i = 3:length(t)-1
    F{i-1} = F_hat{i-1};
    F{i} = F{i-1}(:); F{i} = F{i}';
    F{i} = bsxfun(@max,F{i},S(:,i+1)); F{i} = sort(F{i},2);
    F_hat{i}(:,1) = F{i}(:,1); F_hat{i}(:,end) = F{i}(:,end);
    
    for j = 2:m-1
        F_hat{i}(:,j) = F_hat{i}(:,1) + (F_hat{i}(:,end) - ...
            F_hat{i}(:,1))*(j-1)/(m-1);
    end
    F{i} = F_hat{i};
end

clear F

% LOOKBACK IMMEDIATE PAYOFFS
LB_put = cell(size(F_hat)); V_LB_put = LB_put;

for i = length(LB_put):-1:1
    LB_put{i} = bsxfun(@minus,F_hat{i},S(:,i+1));
end

% LOOKBACK CV
P_augm = zeros(m,m^2);  P_augm(:,1:m) = P{1}; 
for j = 2:m   
    P_augm(j,:) = circshift(P_augm(j,:),[0,m*(j-1)]);
end

V_LB_put{end} = LB_put{end};

for i = length(V_LB_put)-1:-1:2
    V_LB_put{i} = max(LB_put{i},exp(-r*(t(i+2)-t(i+1)))*P{i}*LB_put{i+1});
end

V_LB_put{1} = max(LB_put{1},exp(-r*(t(3)-t(2)))*P_augm*V_LB_put{2}(:));

lookback_put = exp(-r*(t(2)-t(1)))*q'*V_LB_put{1};

% LB CALL

f0 = S(:,1);
f = cell(1,length(t)-1);
f_hat = f; [f_hat{:}] = deal(NaN(m,m));

f{1} = min(f0,S(:,2)); f_hat{1} = f{1};
f{1} = f{1}';
f{2} = bsxfun(@min,f{1},S(:,3)); f_hat{2} = f{2};

for i = 3:length(t)-1
    f{i-1} = f_hat{i-1};
    f{i} = f{i-1}(:); f{i} = f{i}';
    f{i} = bsxfun(@min,f{i},S(:,i+1)); f{i} = sort(f{i},2);
    f_hat{i}(:,1) = f{i}(:,1); f_hat{i}(:,end) = f{i}(:,end);
    
    for j = 2:m-1
        f_hat{i}(:,j) = f_hat{i}(:,1) + (f_hat{i}(:,end) - f_hat{i}(:,1))...
            *(j-1)/(m-1);
    end
    f{i} = f_hat{i};
end

clear f

% LOOKBACK IMMEDIATE PAYOFFS
LB_call = cell(size(f_hat)); V_LB_call = LB_call;

for i = length(LB_call):-1:1
    LB_call{i} = bsxfun(@minus,S(:,i+1),f_hat{i});
end


V_LB_call{end} = LB_call{end};

for i = length(V_LB_call)-1:-1:2
    V_LB_call{i} = max(LB_call{i},exp(-r*(t(i+2)-t(i+1)))*P{i}*LB_call{i+1});
end

V_LB_call{1} = max(LB_call{1},exp(-r*(t(3)-t(2)))*P_augm*V_LB_call{2}(:));

lookback_call = exp(-r*(t(2)-t(1)))*q'*V_LB_call{1};

end

