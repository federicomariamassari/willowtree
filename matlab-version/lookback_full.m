function [lookback_call_example,lookback_put_example] = lookback_full(m,t,S,q,P,r)

% PRODUCES THESIS EXAMPLE (pp.53-58) IF:

% NUMBER OF STATES (m): 3
% LENGTH OF TIME WINDOW (T): 1
% EQUAL TIME STEPS (DIGIT 0)
% NUMBER OF TIME STEPS: 3
% GAMMA = 0.1
% CURRENT PRICE OF THE UNDERLYING (S0) = 10
% LEVEL OF INTEREST RATE (r) = 0.05
% VOLATILITY OF THE UNDERLYING (SIGMA) = 0.3
% q, P, and S are computed according to OPT strategy.

%  [a] CREATE PATH VARIABLE F

F0 = S(:,1);
F = cell(1,length(t)-1); % collects path variables from F{1} to F{end}

F{1} = max(F0,S(:,2));
F{2} = bsxfun(@max,F{1}',S(:,3));

for i = 3:length(t)-1
    F{i} = F{i-1}(:); F{i} = F{i}';
    F{i} = bsxfun(@max,F{i},S(:,i+1));
end

% [b] CREATE EXPANDED P MATRICES (p. 56)

P_expanded = cell(size(P)); [P_expanded{:}] = deal(zeros(m,m^2));
for i = 1:length(P)
    P_expanded{i}(:,1:m) = P{i};

for j = 2:m   
    P_expanded{i}(j,:) = circshift(P_expanded{i}(j,:),[0,m*(j-1)]);
end
end

% [c] MATRICES OF IMMEDIATE PAYOFFS (LB_put)

LB_put = cell(size(F));

for i = length(LB_put):-1:1
    LB_put{i} = bsxfun(@minus,F{i},S(:,i+1));
end

% [d] MATRICES OF OPTION VALUES: COMPARISON BETWEEN IMMEDIATE PAYOFFS AND
%     CONTINUATION VALUES.

V_LB_put = cell(size(F));
V_LB_put{end} = LB_put{end};
V_LB_put{end} = reshape(V_LB_put{end},m^2,m);

for i = length(V_LB_put)-1:-1:2
    V_LB_put{i} = max(LB_put{i},exp(-r*(t(i+2)-t(i+1)))*P_expanded{i}*...
        V_LB_put{i+1});
end

V_LB_put{1} = max(LB_put{1}, exp(-r*(t(3)-t(2)))*P_expanded{1}*...
    V_LB_put{2}(:));

% [e] OPTION PRICE AT VALUATION DATE, t = 0.

lookback_put_example = exp(-r*(t(2)-t(1)))*q'*V_LB_put{1};

%% [2] AMERICAN LOOKBACK CALL OPTION

%  [a] CREATE PATH VARIABLE f

f0 = S(:,1);
f = cell(1,length(t)-1);

f{1} = min(f0,S(:,2));
f{2} = bsxfun(@min,f{1}',S(:,3));

for i = 3:length(t)-1
    f{i} = f{i-1}(:); f{i} = f{i}';
    f{i} = bsxfun(@min,f{i},S(:,i+1));
end

% [c] MATRICES OF IMMEDIATE PAYOFFS (LB_call)

LB_call = cell(size(F));

for i = length(LB_call):-1:1
    LB_call{i} = bsxfun(@minus,S(:,i+1),f{i});
end

% [d] MATRICES OF OPTION VALUES: COMPARISON BETWEEN IMMEDIATE PAYOFFS AND
%     CONTINUATION VALUES.

V_LB_call = cell(size(F));
V_LB_call{end} = LB_call{end};
V_LB_call{end} = reshape(V_LB_call{end},m^2,m);

for i = length(V_LB_call)-1:-1:2
    V_LB_call{i} = max(LB_call{i},exp(-r*(t(i+2)-t(i+1)))*P_expanded{i}*...
        V_LB_call{i+1});
end

V_LB_call{1} = max(LB_call{1}, exp(-r*(t(3)-t(2)))*P_expanded{1}*...
    V_LB_call{2}(:));

% [e] OPTION PRICE AT VALUATION DATE, t = 0.

lookback_call_example = exp(-r*(t(2)-t(1)))*q'*V_LB_call{1};

end

