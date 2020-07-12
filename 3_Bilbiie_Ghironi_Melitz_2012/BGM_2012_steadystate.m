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

Z = 1 ;

% Symmetric equilibrium
Eqs = @(x)...
[...

% (1) Pricing
%p - mu * w/exp(Z);
x(1) - x(2) * x(3)/Z;

% (2-3) Markup & Variety effect
% Translog
   %mu - 1 - 1/(sigm*N);
   %p - exp( -0.5*(N_bar-N)/(sigm*N_bar*N) );
   %x(2) - 1 - 1/(sigm*x(5));
   %x(1) - exp( -0.5*(N_bar-x(5))/(sigm*N_bar*x(5)) );
% CES
   %mu - theta/(theta-1);
   %p - N^(1/(theta-1));
   x(2) - theta/(theta-1);
   x(1) - x(5)^(1/(theta-1));

% (4) Profits
%d = (1-1/mu) * c/N(-1);
x(8) - (1-1/x(2)) * x(4)/x(5);

% (5) Free entry
%v = w * f_e/exp(Z);
x(7) - x(3) * f_e/Z;

% (6) Intratemporal optimality
%chi*L^(1/phi) = w/c;
chi*x(6)^(1/phi) - x(3)/x(4);

% (7) Aggregate accounting
%c + N_e*v = w*L + N(-1)*d;
x(4) + delta/(1-delta)*x(5)*x(7) - x(3)*x(6) - x(5)*x(8);

% (8) Euler equation (shares)
%v = betta * (1-delta) * (v + d);
x(7) - betta * (1-delta) * (x(7) + x(8));

...
];




%-------------------------------------------------------------------------
% Run the solver and save the result
%-------------------------------------------------------------------------

% Initial guesses
%x0 = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];
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

p  = ss(1)
mu = ss(2)
w  = ss(3)
c  = ss(4)
N  = ss(5)
L  = ss(6)
v  = ss(7)
d  = ss(8)

N_e = delta*N/(1-delta);
y = c/(N*p);
Y = c + N_e*v;
D = d*N;
R_E = (v+d)/v;
I_E = N_e*v
L_E = N_e*f_e/Z;
L_C = L - L_E;


log_c = log(c);
log_p = log(p);
log_d = log(d);
log_v = log(v);
log_L = log(L);
log_w = log(w);
log_N = log(N);
log_Z = log(Z);
log_N_e = log(N_e);
log_mu = log(mu);

log_y = log(y);
log_Y = log(Y);
log_D = log(D);
log_R_E = log(R_E);
log_I_E = log(I_E);
log_L_E = log(L_E);
log_L_C = log(L_C);



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
