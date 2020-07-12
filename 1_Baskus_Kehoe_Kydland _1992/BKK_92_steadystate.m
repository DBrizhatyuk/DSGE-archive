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

q_H = 1;
q_F = 1;
A = 1;
A_H = 1;
A_F = 1;
nx_H = 0;
nx_F = 0;
sdf_H = betta;
sdf_F = betta;

% Symmetric equilibrium
Eqs = @(x)...
[...

% Resource constraint
%c + delta*k  - A * k^(theta)*n^(1-theta);
x(1) + delta*x(2)  - exp(A) * x(2)^(theta)*x(3)^(1-theta);

% Labor
%((1-muu)/muu)*(c/(1-n)) - (1-theta)*A * (k/n)^theta;
((1-muu)/muu)*(x(1)/(1-x(3))) - (1-theta)*exp(A) * (x(2)/x(3))^theta;

% Capital
%betta * ( (1-delta) +  theta*A * (k/n)^(theta-1) )  - 1;
betta * ( (1-delta) +  theta*exp(A) * (x(2)/x(3))^(theta-1) )  - 1;

...
];




%-------------------------------------------------------------------------
% Run the solver and save the result
%-------------------------------------------------------------------------

% Initial guesses
%x0 = [1 5 0.5];
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

c_H = ss(1);
c_F = c_H;
k_H = ss(2);
k_F = k_H;
n_H = ss(3);
n_F = n_H;

l_H = 1-n_H;
l_F = l_H;

I_H = delta*k_H;
I_F = I_H;

y_H = c_H + I_H;
y_F = y_H;

w_H = (1-theta)*exp(A_H) * (k_H/n_H)^theta;
w_F = w_H;
R_k_H = theta*exp(A_H) * (k_H/n_H)^(theta-1);
R_k_F = R_k_H;

MU_c_H = muu*(c_H^B1)*(l_H^B2);
MU_c_F = MU_c_H;

log_y_H = log(y_H);
log_y_F = log(y_F);
log_c_H = log(c_H);
log_c_F = log(c_F);
log_I_H = log(I_H);
log_I_F = log(I_F);
log_k_H = log(k_H);
log_k_F = log(k_F);
log_n_H = log(n_H);
log_n_F = log(n_F);
log_A_H = log(A_H);
log_A_F = log(A_F);



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
