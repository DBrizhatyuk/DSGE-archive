function [ys,check] = model_steadystate(ys,exo)

% function [ys,check] = *model*_steadystate(ys,exo)
% computes the steady state for the associated *model*.mod file using a numerical solver
%
% Inputs:
%   - ys        [vector] vector of initial values for the steady state of the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output:
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                        1 of not (allows to impose restriction on parameters)

global M_  % Dynare's model setup
check = 0; % set the indicator vaiable (=1 means that the steady state is not found)




%-------------------------------------------------------------------------
% Upload initial values of structural parameters
%-------------------------------------------------------------------------
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end




%-------------------------------------------------------------------------
% Set up the problem (try to use closed-form solutions as much as possible)
%-------------------------------------------------------------------------

Y = Y_ss; % steady-state TFP calibration
L = 1; % weight of leisure in utility calibration

markup  = epsilon/(epsilon-1);
P       = 1/markup;
R       = 1/betta;
delta_K = delta_K_ss;

% Symmetric equilibrium
Eqs = @(x)...
[...

% (2) Labor supply
% W = chi*L^(1/psi);
% W = chi*L^(1/psi) * (1-h*betta)^(-1)*(C*(1-h))^(gama);
%x(2) - x(3)*L^(1/psi);
x(2) - x(3)*L^(1/psi) * (1-h*betta)^(-1)*(x(1)*(1-h))^(gama);

% (3) Labor demand
% W = P*(1-alphaa)*Y/L;
x(2) - P*(1-alphaa)*Y/L;

% (3) Capital demand
% R_K = (P*alphaa*Y/K + 1 - o.delta_K);
x(4) - (P*alphaa*Y/x(10) + 1 - delta_K);

% (15) Output
% Y = A * K^alphaa * L^(1-alphaa);
Y - x(11) * x(10)^alphaa * L^(1-alphaa);

% (8) Spread
% betta * Omega * (R_K - R) = theta * mu;
betta * x(6) * (x(4) - R) - theta * x(5);

% (9) Ex ante MV of net worth
%  Omega = (1-delta_B) + delta_B*nu;
x(6) - (1-delta_B) - delta_B*x(7);

% (10) Ex post MV of net worth
% nu = (betta*Omega*R)/(1-mu);
x(7) - (betta*x(6)*R)/(1-x(5));

% (11) Net worth
% N = delta_B * ( (R_K-R)*S + R*N ) + omega*S;
x(8) - delta_B * ( (x(4)-R)*x(9) + R*x(8) ) - omega*x(9);

% (12) Binding constraint
% nu*N = theta*S;
x(7)*x(8) - theta*x(9);

% (13) Banks resourse constraint
% B + N = S;
x(12) + x(8) - x(9);

% (14) Final output resourse constraint
% Y = C + delta_K*K;
Y - x(1) - delta_K*x(10);

% (15)
% K = S;
x(10) - x(9);

...
];




%-------------------------------------------------------------------------
% Run the solver and save the result
%-------------------------------------------------------------------------

% Initial guesses
%x0 =  [ 0.8006    0.6091    0.6091    1.0126    0.0101    1.5457    1.5614   1.946    7.9771    7.9771    0.5040    6.0306];
x0 = load('steady_st_init_values');

% Fsolve options
tol      = 1e-10;
Iter     = 100000;
FunEvals = 100000;
options  = optimset('Display', 'on', 'MaxFunEvals', FunEvals, 'MaxIter', Iter, 'TolFun', tol, 'TolX', tol);

% Run fsolve
[ss,fval,exitflag] = fsolve(Eqs, x0.ss, options);

% Exitflag indicates the reason fsolve stopped;  exitflag < 1 means that the solution was not found
if exitflag <1
    check=1;
    return;
end

% Save the solution
save('steady_st_init_values.mat','ss');




%-------------------------------------------------------------------------
% Assign output to model variables
%-------------------------------------------------------------------------

ss = real(ss);

C         = ss(1);       % consumption
W         = ss(2);       % real wage
chi       = ss(3);       % labor disutility
R_K       = ss(4);       % return on capital
mu        = ss(5);       %
Omega     = ss(6);       %
nu        = ss(7);       %
N         = ss(8);       % aggregate banks' net worth
S         = ss(9);       % aggregate deposits
K         = ss(10);      % capital
A         = ss(11);      % steady-state productivity
B         = ss(12);      % deposits

lambda = (1-h*betta) * (C*(1-h))^(-gama); % Separable utility
%lambda = ( C - h*C - chi*L^(1+1/psi)/(1+1/psi) )^(-gama) - h*( C - h*C - chi*L^(1+1/psi)/(1+1/psi) )^(-gama); % GHH utility

I = delta_K*K;
I_ss = I;
I_n = I - delta_K*K;
q = 1;
U = 1;
c1 = P*alphaa*Y/K;

i = i_ss;
P = (epsilon-1)/epsilon;
pi = pi_ss;
sdf = betta;
sdf_adjusted = betta*Omega;
spread = R_K - R;
xi = 1;

leverage = S*q/N;

log_xi = log(xi);
log_Y = log(Y);
log_C = log(C);
log_I = log(I);
log_K = log(K);
log_L = log(L);
log_q = log(q);
log_N = log(N);

%-------------------------------------------------------------------------
% Save the structural parameters (the calibrated ones will be updated)
%-------------------------------------------------------------------------

for iter = 1:length(M_.params) % update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ]);
end




%-------------------------------------------------------------------------
% Save the calibrated steady state
%-------------------------------------------------------------------------

% constructs the output vestor ys
NumberOfEndogenousVariables = M_.orig_endo_nbr; % auxiliary variables are set automatically

steady_st_valuses = zeros(1,NumberOfEndogenousVariables);

for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
  eval([ 'steady_st_values(' num2str(ii) ') = ' varname ';']);
end

save('steady_st_values.mat','steady_st_values');
