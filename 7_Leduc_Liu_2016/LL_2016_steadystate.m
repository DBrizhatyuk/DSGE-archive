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

Z = 1;
sigma_Z = sigma_Z_ss;
i = i_ss;
p = (eta-1)/eta;

Eqs = @(x)...
[...

% (1) Market clearing
%Z*L = c + kappa*v;
Z*x(1) - x(2) - kappa*x(3);

% (2) Law of motion for labor
%sigm*L = q*v;
sigm*x(1) - x(4)*x(3);

% (3) Job creation
%(1-betta*(1-sigm)) * kappa/q = p*Z - w;
(1-betta*(1-sigm)) * kappa/x(4) - p*Z + x(5);

% (4) Wage curve w/ real wage rigidity
%w = (1-phi)*(chi/( (1-betta*h)*(c - h*c)^(-gama) ) + u) + phi*(p*Z + betta*(1-sigm)*theta*kappa);
x(5) - (1-phi)*(chi/( (1-betta*h)*(x(2) - h*x(2))^(-gama) ) + u) - phi*(p*Z + betta*(1-sigm)*x(6)*kappa);

% (5) Vacancy filling rate
%q = chi_M*theta^(-epsilon);
x(4) - chi_M*x(6)^(-epsilon);

% (6) [Inverse of] labor market tightness
%theta = v/(1-(1-sigm)*L);
x(6) - x(3)/(1-(1-sigm)*x(1));

...
];



%-------------------------------------------------------------------------
% Run the solver and save the result
%-------------------------------------------------------------------------

% Initial guesses
%x0 = [0.5 0.5 0.5 0.5 0.5 0.5];
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

L = ss(1)
c = ss(2)
v = ss(3)
q = ss(4)
w = ss(5)
theta = ss(6)

Y = Z*L;
Y_ss = Y;
lambda = (1-betta*h)*(c - h*c)^(-gama);
w_N = w;
U = 1 - (1-sigm)*L - q*v;
pi = pi_ss;
sdf = betta;
R = (1+i)/pi;
J = kappa/q;

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
