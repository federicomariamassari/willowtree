fprintf(' \n <strong>A WILLOW TREE LATTICE TO PRICE PATH-DEPENDENT OPTIONS</strong>\n')
fprintf(' Massari Federico Maria\n')

m = input('\n Number of states (m)? ');

T = input('\n Length of time window (in years)? ');

fprintf(' \n Choice of time steps. Digit 0 for EQUAL, or 1 for UNEQUAL.\n'); 
fprintf(' \n If UNEQUAL, please provide row vector t. Either call a vector from\n');
fprintf(' the workspace, and digit the name of the variable, or input manually.\n')
fprintf(' In the latter case, digit vector without the name of the variable.\n')
fprintf(' Example: digit [t(1) t(2) ... t(end)] and NOT t = [t(1) t(2) ... t(end)].\n')
fprintf(' To obtain meaningful output, t should be bound in [0,T], so t(1) = 0,\n');
fprintf(' and t(end) = T. Example: T = 0.5, t = [0 0.11 0.28 0.41 0.5].\n');
time_steps = input('\n EQUAL or UNEQUAL time steps? (0=EQUAL, 1=UNEQUAL)? ');

if time_steps == 0
    t = input('\n Number of time steps? ');
    t = 0:T/t:T;
elseif time_steps == 1
    t = input('\n Manually enter the (row) vector of time steps: ');
end

% [A] SAMPLING STRATEGIES
fprintf('\n\n <strong>*********************************************************************</strong>');
fprintf('\n <strong>[A] SAMPLING STRATEGIES</strong>\n');
fprintf(' <strong>*********************************************************************</strong>\n');

% CURRAN'S STRATEGY

fprintf('\n <strong>[1] CURRAN''S METHODOLOGY: {z(i),q(i)}</strong>\n');
[z_CUR, q_CUR, G_CUR] = sampling(m,t,0,'Curran');

fprintf('\n <strong>STATISTICS</strong>\n');
fprintf('\n Value of mean = %g\n', q_CUR'*z_CUR);
fprintf(' Value of variance = %g\n', q_CUR'*z_CUR.^2);
fprintf(' Value of kurtosis = %g\n', kurtosis(z_CUR));
fprintf('\n  (a) When the pairs {z(i),q(i)} are selected according to Curran''s\n');
fprintf('      methodology, the LP problem has no solution. The vectors cannot\n');
fprintf('      be chosen arbitrarily, but in such a way to create some linear\n');
fprintf('      dependence (LD) between the matrix of constraints (Aeq) and the\n');
fprintf('      vector of solutions (beq). If not, rank([Aeq beq]) > rank(Aeq),\n');
fprintf('      and the system has no unique solution.\n');
fprintf('\n  (b) When m is small, the choice of {z(i),q(i)} produces values of\n');
fprintf('      variance and kurtosis not in line with those of the standard\n');
fprintf('      normal distribution.\n');

% IMPROVED CURRAN'S METHODOLOGY (OPT STRATEGY WITH GAMMA = 0)

fprintf('\n <strong>[2] IMPROVED CURRAN''S METHODOLOGY: {z(i),q(i)}</strong>\n');
fprintf('\n This methodology selects the vector of standard normal variates <strong>z</strong>\n');
fprintf(' to make kurtosis closer to 3. It is OPT strategy with gamma = 0.\n');

fprintf('\n <strong>CONSTRAINED MINIMIZATION PROBLEM</strong>\n');
[z_iCUR, q_iCUR, G_iCUR] = sampling(m,t,0,'new sampling strategy','OPT strategy');

fprintf('\n <strong>STATISTICS</strong>\n');
fprintf('\n Value of mean = %g\n', q_iCUR'*z_iCUR);
fprintf(' Value of variance = %g\n', q_iCUR'*z_iCUR.^2);
fprintf(' Value of kurtosis* = %g\n', kurtosis(z_iCUR));
fprintf('\n * When gamma is equal to zero, kurtosis is calculated according\n');
fprintf('   to the built-in MATLAB formula ''kurtosis''.\n');

fprintf('\n <strong>LINEAR PROGRAMMING PROBLEM</strong>\n');
[P_iCUR,fval_c,exitflag_c,h_iCUR] = LP(m,z_iCUR,q_iCUR,t,0,'new sampling strategy');

% OPTIMISATION STRATEGY

fprintf('\n <strong>[3a] NEW SAMPLING STRATEGY: OPT STRATEGY FOR {z(i),q(i)}</strong>\n');

gamma = input('\n Value of gamma? ');

fprintf('\n <strong>CONSTRAINED MINIMIZATION PROBLEM</strong>\n');
[z_OPT, q_OPT, G_OPT] = sampling(m,t,gamma,'new sampling strategy','OPT strategy');

fprintf('\n <strong>STATISTICS</strong>\n');
fprintf('\n Value of mean = %g\n', q_OPT'*z_OPT);
fprintf(' Value of variance = %g\n', q_OPT'*z_OPT.^2);
fprintf(' Value of kurtosis* = %g\n', q_OPT'*z_OPT.^4);
fprintf('\n * When gamma is different from zero, kurtosis is calculated according\n');
fprintf('   to the following formula: q''*z.^4\n');
fprintf('   Otherwise, according to the built-in MATLAB formula ''kurtosis''.\n');
fprintf('\n <strong>LINEAR PROGRAMMING PROBLEM</strong>\n');
[P_OPT,fval_n,exitflag_n,h_OPT] = LP(m,z_OPT,q_OPT,t,gamma,'new sampling strategy');

% ADJUSTED PROBABILITY STRATEGY

if mod(m,2) == 0
fprintf('\n <strong>[3b] NEW SAMPLING STRATEGY: ADJ STRATEGY FOR {z(i),q(i)}</strong>\n');

[z_ADJ, q_ADJ, G_ADJ] = sampling(m,t,gamma,'new sampling strategy','ADJ strategy');

fprintf('\n <strong>STATISTICS</strong>\n');
fprintf('\n Value of mean = %g\n', q_ADJ'*z_ADJ);
fprintf(' Value of variance = %g\n', q_ADJ'*z_ADJ.^2);
fprintf(' Value of kurtosis* = %g\n', q_ADJ'*z_ADJ.^4);
fprintf('\n * When gamma is different from zero, kurtosis is calculated according\n');
fprintf('   to the following formula: q''*z.^4\n');
fprintf('   Otherwise, according to the built-in MATLAB formula ''kurtosis''.\n');
fprintf('\n  (a) When m is small, the choice of {z(i),q(i)} produces a value of\n');
fprintf('      variance not in line with that of the standard normal distribution;\n');
fprintf('      the value of kurtosis is consistent due to the optimization procedure.\n');
end

% FIRST PARTIAL MOMENT MATCHING STRATEGY

fprintf('\n <strong>[3b] NEW SAMPLING STRATEGY: FPM STRATEGY FOR {z(i),q(i)}</strong>\n');

fprintf('\n <strong>CONSTRAINED MINIMIZATION PROBLEM</strong>\n');
[z_FPM, q_FPM, G_FPM] = sampling(m,t,gamma,'new sampling strategy','FPM strategy');

fprintf('\n <strong>STATISTICS</strong>\n');
fprintf('\n Value of mean = %g\n', q_FPM'*z_FPM);
fprintf(' Value of variance = %g\n', q_FPM'*z_FPM.^2);
fprintf(' Value of kurtosis* = %g\n', q_FPM'*z_FPM.^4);
fprintf('\n * When gamma is different from zero, kurtosis is calculated according\n');
fprintf('   to the following formula: q''*z.^4.\n');
fprintf('   Otherwise, according to the built-in MATLAB formula ''kurtosis''.\n');
fprintf('\n  (a) When m is small, the choice of {z(i),q(i)} produces a value of\n');
fprintf('      kurtosis not in line with that of the standard normal distribution;\n');
fprintf('      the value of variance is consistent due to the optimization procedure.\n');

fprintf('\n <strong>LINEAR PROGRAMMING PROBLEM</strong>\n');
[P_FPM,fval_FPM,exitflag_FPM,h_FPM] = LP(m,z_FPM,q_FPM,t,gamma,'new sampling strategy');

% [B] PRICING OF EUROPEAN, PLAIN VANILLA DERIVATIVES
fprintf('\n\n <strong>*********************************************************************</strong>');
fprintf('\n <strong>[B] PRICING OF EUROPEAN, PLAIN VANILLA DERIVATIVES</strong>\n');
fprintf(' <strong>*********************************************************************</strong>\n');

S0 = input('\n Current price of the underlying (S0)? ');
K = input('\n Strike price (K)? ');
r = input('\n Level of interest rate (r)? ');
sigma = input('\n Volatility of the underlying (sigma)? ');

fprintf('\n <strong>[0] BLACK-SCHOLES-MERTON CLOSED FORMULA</strong>\n');
[call_BSM, put_BSM] = blsprice(S0,K,r,T,sigma);
fprintf('\n True price of European call option = %g\n', call_BSM);
fprintf(' True price of European put option = %g\n', put_BSM);

fprintf('\n <strong>[1] CURRAN''S METHODOLOGY</strong>\n');
fprintf(' ');
[S_CUR, call_CUR, put_CUR] = pricing(S0,K,r,sigma,m,q_CUR,G_CUR,[],t,T,'CUR');
[error_call_CUR,error_put_CUR] = pricing_error(call_BSM,call_CUR,put_BSM,put_CUR);
fprintf('\n Value of European call option = %g\n', call_CUR);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_CUR);
fprintf('\n Value of European put option = %g\n', put_CUR);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_CUR);

fprintf('\n <strong>[2] IMPROVED CURRAN''S METHODOLOGY</strong>\n');
fprintf(' ');
[S_iCUR, call_iCUR, put_iCUR] = pricing(S0,K,r,sigma,m,q_iCUR,G_iCUR,P_iCUR,t,T,'iCUR');
[error_call_iCUR,error_put_iCUR] = pricing_error(call_BSM,call_iCUR,put_BSM,put_iCUR);
fprintf('\n Value of European call option = %g\n', call_iCUR);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_iCUR);
fprintf('\n Value of European put option = %g\n', put_iCUR);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_iCUR);

fprintf('\n <strong>[3a] NEW SAMPLING STRATEGY: OPT STRATEGY</strong>\n');
fprintf(' ');
[S_OPT, call_OPT, put_OPT] = pricing(S0,K,r,sigma,m,q_OPT,G_OPT,P_OPT,t,T,'OPT');
[error_call_OPT,error_put_OPT] = pricing_error(call_BSM,call_OPT,put_BSM,put_OPT);
fprintf('\n Value of European call option = %g\n', call_OPT);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_OPT);
fprintf('\n Value of European put option = %g\n', put_OPT);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_OPT);

if mod(m,2) == 0
fprintf('\n <strong>[3b] NEW SAMPLING STRATEGY: ADJ STRATEGY</strong>\n');
fprintf(' ');
[S_ADJ, call_ADJ, put_ADJ] = pricing(S0,K,r,sigma,m,q_ADJ,G_ADJ,[],t,T,'ADJ');
[error_call_ADJ,error_put_ADJ] = pricing_error(call_BSM,call_ADJ,put_BSM,put_ADJ);
fprintf('\n Value of European call option = %g\n', call_ADJ);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_ADJ);
fprintf('\n Value of European put option = %g\n', put_ADJ);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_ADJ);
end

fprintf('\n <strong>[3c] NEW SAMPLING STRATEGY: FPM STRATEGY</strong>\n');
fprintf(' ');
[S_FPM, call_FPM, put_FPM] = pricing(S0,K,r,sigma,m,q_FPM,G_FPM,P_FPM,t,T,'FPM');
[error_call_FPM,error_put_FPM] = pricing_error(call_BSM,call_FPM,put_BSM,put_FPM);
fprintf('\n Value of European call option = %g\n', call_FPM);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_FPM);
fprintf('\n Value of European put option = %g\n', put_FPM);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_FPM);

fprintf('\n <strong>[4] BINOMIAL PRICING MODEL WITH m STEPS</strong>\n');
fprintf('\n NOTE: The code is a modified version of listing "euro5.m", in\n');
fprintf(' Higham, D.J.: "Nine Ways To Implement The Binomial Method for Option\n');
fprintf(' Valuation in MATLAB", SIAM Review, 44(4), 2002, 661-677.\n');
fprintf('\n ');
[call_BIN, put_BIN] = binomial(S0,K,r,sigma,T,t,'m');
[error_call_BIN,error_put_BIN] = pricing_error(call_BSM,call_BIN,put_BSM,put_BIN);
fprintf('\n Value of European call option = %g\n', call_BIN);
fprintf(' Relative Error to BSM call formula = %g\n', error_call_BIN);
fprintf('\n Value of European put option = %g\n', put_BIN);
fprintf(' Relative Error to BSM put formula = %g\n', error_put_BIN);

% [C] PRICING OF AMERICAN, PLAIN VANILLA DERIVATIVES

fprintf('\n\n <strong>*********************************************************************</strong>');
fprintf('\n <strong>[C] PRICING OF AMERICAN, PLAIN VANILLA DERIVATIVES</strong>\n');
fprintf(' <strong>*********************************************************************</strong>\n');

fprintf('\n <strong>[0] BINOMIAL PRICING MODEL (5000 STEPS)</strong>\n');
fprintf('\n The binomial model with 5000 steps is chosen as a reliable proxy\n');
fprintf(' for the true price of American, plain vanilla call and put options.\n');
fprintf(' The relative errors of the willow tree strategies are computed using\n');
fprintf(' the binomial model figures as benchmark.\n');
fprintf('\n NOTE: The code is a modified version of listing "american.m", in\n');
fprintf(' Higham, D.J.: "Nine Ways To Implement The Binomial Method for Option\n');
fprintf(' Valuation in MATLAB", SIAM Review, 44(4), 2002, 661-677.\n');

[call_BIN_American, put_BIN_American] = binomial_american(S0,K,r,sigma,T);
fprintf('\n Value of American call option = %g\n', call_BIN_American);
fprintf(' Value of American put option = %g\n', put_BIN_American);
fprintf(' ');

[call_BIN_5000, put_BIN_5000] = binomial(S0,K,r,sigma,T,t,'5000');
fprintf('\n Value of European call option = %g\n', call_BIN_5000);
fprintf(' Value of European put option = %g\n', put_BIN_5000);


fprintf('\n <strong>[2] IMPROVED CURRAN''S METHODOLOGY</strong>\n');
fprintf(' ');
[call_iCUR_American, put_iCUR_American] = pricing_american(m,q_iCUR,t,G_iCUR,P_iCUR,S0,K,r,sigma,'iCUR');
[error_call_iCUR_American,error_put_iCUR_American] = pricing_error_american(call_BIN_American,call_iCUR_American,put_BIN_American,put_iCUR_American);
fprintf('\n Value of American call option = %g\n', call_iCUR_American);
fprintf(' Relative Error to BIN call formula = %g\n', error_call_iCUR_American);
fprintf('\n Value of American put option = %g\n', put_iCUR_American);
fprintf(' Relative Error to BIN put formula = %g\n', error_put_iCUR_American);

fprintf('\n <strong>[3a] NEW SAMPLING STRATEGY: OPT STRATEGY</strong>\n');
fprintf(' ');
[call_OPT_American, put_OPT_American] = pricing_american(m,q_OPT,t,G_OPT,P_OPT,S0,K,r,sigma,'OPT');
[error_call_OPT_American,error_put_OPT_American] = pricing_error_american(call_BIN_American,call_OPT_American,put_BIN_American,put_OPT_American);
fprintf('\n Value of American call option = %g\n', call_OPT_American);
fprintf(' Relative Error to BIN call formula = %g\n', error_call_OPT_American);
fprintf('\n Value of American put option = %g\n', put_OPT_American);
fprintf(' Relative Error to BIN put formula = %g\n', error_put_OPT_American);

fprintf('\n <strong>[3c] NEW SAMPLING STRATEGY: FPM STRATEGY</strong>\n');
fprintf(' ');
[call_FPM_American, put_FPM_American] = pricing_american(m,q_FPM,t,G_FPM,P_FPM,S0,K,r,sigma,'FPM');
[error_call_FPM_American,error_put_FPM_American] = pricing_error_american(call_BIN_American,call_FPM_American,put_BIN_American,put_FPM_American);
fprintf('\n Value of American call option = %g\n', call_FPM_American);
fprintf(' Relative Error to BIN call formula = %g\n', error_call_FPM_American);
fprintf('\n Value of American put option = %g\n', put_FPM_American);
fprintf(' Relative Error to BIN put formula = %g\n', error_put_FPM_American);

% [D] PRICING OF AMERICAN, LOOKBACK DERIVATIVES

fprintf('\n\n <strong>*********************************************************************</strong>');
fprintf('\n <strong>[D] PRICING OF AMERICAN, LOOKBACK DERIVATIVES</strong>\n');
fprintf(' <strong>*********************************************************************</strong>\n');

fprintf('\n <strong>EXPLANATORY PURPOSE ONLY! VERY INACCURATE RESULTS.</strong>\n');
fprintf('\n The results on path dependent options are for explanatory purpose only,\n');
fprintf(' and are not supposed to be accurate, due to the introduction of a set\n');
fprintf(' of auxiliary vectors of representative running maxima, in the forward\n');
fprintf(' shooting grid methodology.\n');
fprintf('\n NOTE: If the chosen number of time steps is t = 3, the code also allows\n');
fprintf(' to compute the prices of lookback options without recourse to auxiliary\n');
fprintf(' vectors of representative running maxima ("full lookback prices").\n');
fprintf(' Due to the very low number of time steps, these results are not supposed\n');
fprintf(' to be accurate either.\n');

if numel(t)-1 == 3
X = input('\n Compute the full lookback prices? (0=NO, 1=YES) ');
if X == 1
    
    fprintf('\n <strong>[A] FULL LOOKBACK PRICES</strong>\n');
    
    fprintf('\n <strong>[2] IMPROVED CURRAN''S METHODOLOGY</strong>\n');
    fprintf(' ');
    [LB_call_example_iCUR,LB_put_example_iCUR] = lookback_full(m,t,S_iCUR,q_iCUR,P_iCUR,r);
    fprintf('\n Value of American lookback call option = %g\n', LB_call_example_iCUR);
    fprintf(' Value of American lookback put option = %g\n', LB_put_example_iCUR);
    
    fprintf('\n <strong>[3a] NEW SAMPLING STRATEGY: OPT STRATEGY</strong>\n');
    fprintf(' ');
    [LB_call_example_OPT,LB_put_example_OPT] = lookback_full(m,t,S_OPT,q_OPT,P_OPT,r);
    fprintf('\n Value of American lookback call option = %g\n', LB_call_example_OPT);
    fprintf(' Value of American lookback put option = %g\n', LB_put_example_OPT);
    
    fprintf('\n <strong>[3c] NEW SAMPLING STRATEGY: FPM STRATEGY</strong>\n');
    fprintf(' ');
    [LB_call_example_FPM,LB_put_example_FPM] = lookback_full(m,t,S_FPM,q_FPM,P_FPM,r);
    fprintf('\n Value of American lookback call option = %g\n', LB_call_example_FPM);
    fprintf(' Value of American lookback put option = %g\n', LB_put_example_FPM);
    
    fprintf('\n <strong>[B] WITH AUXILIARY VECTORS OF REPRESENTATIVE RUNNING MAXIMA</strong>\n');
end
end

fprintf('\n <strong>[2] IMPROVED CURRAN''S METHODOLOGY</strong>\n');
fprintf(' ');
[LB_put_iCUR, LB_call_iCUR] = lookback_american(m,t,S_iCUR,q_iCUR,P_iCUR,r);
fprintf('\n Value of American lookback call option = %g\n', LB_call_iCUR);
fprintf(' Value of American lookback put option = %g\n', LB_put_iCUR);

fprintf('\n <strong>[3a] NEW SAMPLING STRATEGY: OPT STRATEGY</strong>\n');
fprintf(' ');
[LB_put_OPT, LB_call_OPT] = lookback_american(m,t,S_OPT,q_OPT,P_OPT,r);
fprintf('\n Value of American lookback call option = %g\n', LB_call_OPT);
fprintf(' Value of American lookback put option = %g\n', LB_put_OPT);

fprintf('\n <strong>[3c] NEW SAMPLING STRATEGY: FPM STRATEGY</strong>\n');
fprintf(' ');
[LB_put_FPM, LB_call_FPM] = lookback_american(m,t,S_FPM,q_FPM,P_FPM,r);
fprintf('\n Value of American lookback call option = %g\n', LB_call_FPM);
fprintf(' Value of American lookback put option = %g\n', LB_put_FPM);

% [E] GRAPHICAL REPRESENTATION
fprintf('\n\n <strong>*********************************************************************</strong>');
fprintf(' \n <strong>[E] GRAPHICAL REPRESENTATION OF THE WILLOW TREE</strong>\n')
fprintf(' <strong>*********************************************************************</strong>\n');

fprintf(' \n Works best for a small number of nodes and time steps. Optimised for\n')
fprintf(' m = 11 and t = 4.\n')

X = input('\n Graphical representation? (0=NO, 1=YES) ');
if X == 1
    
fprintf(' \n <strong>FIGURE 1: IMPROVED CURRAN''S METHODOLOGY.</strong> ');
graph_willow(m,z_iCUR,q_iCUR,t,h_iCUR,T,G_iCUR,P_iCUR)
fprintf(' \n <strong>FIGURE 2: OPT STRATEGY.</strong> ');
graph_willow( m,z_OPT,q_OPT,t,h_OPT,T,G_OPT,P_OPT )
fprintf(' \n <strong>FIGURE 3: FPM STRATEGY.</strong> ');
graph_willow( m,z_FPM,q_FPM,t,h_FPM,T,G_FPM,P_FPM )
end

X = input('\n Save output to file? (0=NO, 1=YES) ');
if X == 1
    clear time_steps X
    uisave;
end
clear time_steps X