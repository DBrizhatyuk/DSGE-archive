//--------------------------------------------------------------------------------------------------------------
// Optimal monetary policy under commitment, cost-pust shock in a basic NK model
// Similar to Gali (2008), fig 5.1 & 5.2
//--------------------------------------------------------------------------------------------------------------




var C  // consumption
    Y  log_Y // output
    L  // labor
    mu // markup
    pi // inflation
    i  // nominal interest rate
    u; // markup shock

varexo e_mu;

parameters beta psi_pi theta epsilon alpha
           rho_mu i_ss Y_ss Z;




//-------------------------------------------------------------------------
//               Parameters
//-------------------------------------------------------------------------

beta    = 0.99; // discount factor
psi_pi  = 100;  // price adustment cost
theta   = 6;    // elasticity of substitution
epsilon = 1;    // inverse elasticity of labor supply
alpha   = 1/3;  // production function curvature

i_ss    = 1/beta - 1;
Z       = 1;
Y_ss    = 1;
rho_mu  = 0;

// Competitive economy solution, ad hoc Taylor rule
@#define competitive = 0
//      OR
// Optimal monetary policy under commitment
@#define commitment  = 1




//-------------------------------------------------------------------------
//               Model equilibrium conditions
//-------------------------------------------------------------------------

model;

// (1) Resourse constraint
Y = C + (psi_pi/2)*(pi-1)^2 * Y;

// (2) Output
Y = Z*L^(1-alpha);

// (3) Euler equation
beta*(C/C(+1))* (1+i)/pi(+1) = 1;

// (4) Labor optimality
( (1-alpha)*L^(-alpha)*Z )/mu = L^(epsilon) * C;

// (5) Markup
mu = theta/( (theta-1) + psi_pi * (  pi*(pi-1) - beta*(C/C(+1))*pi(+1)*(pi(+1)-1)*(Y(+1)/Y) )  ) * u;

// (6) Markup shock
log(u) = rho_mu*log(u(-1)) + e_mu;

// (7) Taylor rule
@#if competitive == 1
(1+i)/(1+i_ss) = pi^2;
@#endif

// (8) Log output
log_Y = log(Y);

end;




//-------------------------------------------------------------------------
//               Steady state
//-------------------------------------------------------------------------

steady_state_model;
  mu = theta/(theta-1);
  i = 1/beta - 1;
  L = ((1-alpha)/mu)^(1/(1+epsilon));
  Y = Z*L^(1-alpha);
  C = Y;
  u = 1;
  pi = 1;
  log_Y = log(Y);
end;




//-------------------------------------------------------------------------
//               Simulations
//-------------------------------------------------------------------------

shocks;
var e_mu = 0.01^2;
end;

// IRFs under Taylor rule

@#if competitive == 1
steady;
resid;
check;
stoch_simul(irf=12, order=1, nomoments, nograph, nofunctions);
@#endif

// IRFs under optimal policy under commitment

@#if commitment == 1
planner_objective log(C) - L^(1+epsilon)/(1+epsilon);
ramsey_policy(instruments=(i),irf=12,planner_discount=beta, nograph);
@#endif

// Price level IRF

pi_lvl = 1 + oo_.irfs.pi_e_mu;
price_lvl = ones(1+length(pi_lvl),1);

for i=2:(1+length(pi_lvl))
  price_lvl(i) = price_lvl(i-1) * pi_lvl(i-1);
end
oo_.irfs.price_lvl_e_mu = price_lvl - ones(1+length(pi_lvl),1);

// Plot

varstoplot = {'Y'       'pi'        'price_lvl'   'u'};
varnames   = {'Output'  'Inflation' 'Price level' 'Cost push shock'};

level_vars = {'pi' 'price_lvl'};
irf_len    = 12;
x_tick     = 2;
subplot_x  = 2;
subplot_y  = 2;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);
