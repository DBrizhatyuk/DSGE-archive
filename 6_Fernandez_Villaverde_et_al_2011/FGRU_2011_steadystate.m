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

D = D_ss;
R = R_ss;
e_R = 0;
sigma_R = sigma_R_ss;
q = 1;

Eqs = @(x)...
[...

% (1) Resourse constraint
%A * L^(1-alph) * K^(alph) - C - delta*K = D - D/(1+R);
A * x(2)^(1-alph) * x(1)^(alph) - x(3) - delta*x(1) - D + D/(1+R);

% (2) MU_C
%lambda = ( C - chi*L^(1+psi)/(1+psi) )^(-gamma);
%lambda = C^(-gamma);
%x(4) = ( x(3) - chi*x(2)^(1+psi)/(1+psi) )^(-gamma);
x(4) - x(3)^(-gamma);

% (3) Capital mrkt equilibrium
%1 = beta * ( 1-delta + alph*A*L^(1-alph)*K^(alph-1) );
1 - beta * ( 1-delta + alph*A*x(2)^(1-alph)*x(1)^(alph-1) );

% (4) Labor mrkt equilibrium
%chi*L^(psi) = (1-alph) * A*L^(-alph)*K^(alph);
%chi*L^(psi) / lambda = (1-alph) * A*L^(-alph)*K^(alph);
%hi*x(2)^(psi) = (1-alph) * A*x(2)^(-alph)*x(1)^(alph);
chi*x(2)^(psi) / x(4) - (1-alph) * A*x(2)^(-alph)*x(1)^(alph);
...
];



%-------------------------------------------------------------------------
% Run the solver and save the result
%-------------------------------------------------------------------------

% Initial guesses
%x0 = [0.5 0.5 0.5 0.5];
x0 = load('steady_st_init_values');

% Fsolve options
tol      = 1e-5;
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

K = ss(1)
L = ss(2)
C = ss(3)
lambda = ss(4)

CA = 0;
Y = A * L^(1-alph) * K^(alph);
I = delta*K;
NX = Y - C - I;




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
