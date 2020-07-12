//-----------------------------------------------------------------------------------
// Leduc and Liu (2016): Uncertainty Shocks are Aggregate Demand Shocks
//-----------------------------------------------------------------------------------




var Y // output
    L // aggregate employment
    U // unemployment rate
    c // consumption
    v //
    q //
    theta //
    lambda //
    sdf //
    w   //
    w_N //
    p   //
    pi  //
    i   // nominal interest rate
    R   // real interest rate
    J   // matrch value
    Z sigma_Z; // productivity process

varexo e_Z e_sigma_Z;

parameters sigm betta phi epsilon gama eta omega pi_ss Y_ss i_ss chi h rho_w phi_pi phi_Y kappa chi_M u //
           rho_Z sigma_Z_ss rho_sigma_Z sigma_sigma_Z;    // productivity shock parameters




//-------------------------------------------------------------------------
//               Parameters
//-------------------------------------------------------------------------
gama   = 1;      // elasticity of intertemporal substitution
betta    = 0.99;   // discount rate
chi     = 0.547;  // scale of disutility of working
h       = 0;      // consumption habit
eta     = 10;     // elasticity of substitution b/w intermediate goods
epsilon = 0.5;    // matching function parameter
chi_M   = 0.645;  // matching efficiency
sigm   = 0.1;    // job destruction rate
u       = 0.25;   // flow benefit of unemployment
kappa   = 0.14;   // flow cost of vacancy
phi     = 0.5;    // Nash barganing weight
rho_w   = 0.8;    // real wage rigidity
omega   = 112;    // price adjustment cost

pi_ss   = 1.005;  // 2% annual inflation
i_ss    = pi_ss/betta - 1;
phi_pi  = 1.5;    // Taylor rule, inflation
phi_Y   = 0.2;    // Taylor rule, output
Y_ss    = 1;      // steady-state GDP

// Productivity shock parameter:
rho_Z         = 0.95;
sigma_Z_ss    = 0.01;
rho_sigma_Z   = 0.76;
sigma_sigma_Z = 0.392;




//-------------------------------------------------------------------------
//               Model equilibrium conditions
//-------------------------------------------------------------------------

model;

// (1) Output
Y = Z*L;

// (2) Market clearing
Y*( 1 - (omega/2)*((pi/pi_ss)-1)^2 ) = c + kappa*v;

// (3) Law of motion for labor
L = (1-sigm)*L(-1) + q*v;

// (4) Unemployment rate
U = 1 - (1-sigm)*L(-1) - q*v;

// (5) Job creation
kappa/q = p*Z - w + sdf*(1-sigm)*kappa/q(+1) ;

// (6)-(7) Wage curve w/ real wage rigidity
w_N = (1-phi)*(chi/lambda + u) + phi*(p*Z + sdf*(1-sigm)*theta(+1)*kappa);
w = w(-1)^(rho_w) * w_N^(1-rho_w);

// (8) Vacancy filling rate
q = chi_M*theta^(-epsilon);

// (9) [Inverse of] labor market tightness
theta = v/(1-(1-sigm)*L(-1));

// (10) Relative price of an intermadiate good
p = (eta-1)/eta + (omega/eta) * ( (pi/pi_ss-1)*pi/pi_ss - sdf*(Y(+1)/Y)*(pi(+1)/pi_ss-1)*pi(+1)/pi_ss );

// (11) Consumption Euler
sdf*(1+i)/pi(+1) = 1;

// (12) Taylor rule
(1+i) = (1+i_ss) * (pi/pi_ss)^(phi_pi) *  (Y/Y_ss)^(phi_Y);

// (13) Real interest rate
R = (1+i)/pi(+1);

// (14) Match value
J = kappa/q;

// (15) MU of C
lambda = (c - h*c(-1))^(-gama) - betta*h*(c(+1) - h*c)^(-gama);

// (16) Stochastic discount factor
sdf = betta*(lambda(+1)/lambda);

//-------------------------------------------------------------------------

// (17-18)  Productivity shock
log(Z)       =  rho_Z*log(Z(-1)) + e_Z*sigma_Z;
log(sigma_Z) = (1-rho_sigma_Z)*log(sigma_Z_ss) + rho_sigma_Z*log(sigma_Z(-1)) + e_sigma_Z*sigma_sigma_Z;

end;




//-------------------------------------------------------------------------
//               Steady state
//-------------------------------------------------------------------------

steady;
resid(1);




//-------------------------------------------------------------------------
//               Impulse responses
//-------------------------------------------------------------------------

shocks;
var e_Z = 1;
var e_sigma_Z = 1;
end;


stoch_simul(order=3,periods=0,irf=0,noprint,nograph,nomoments,nofunctions,nocorr,pruning);

irfburninlength   = 5000;
irflength         = 20;

oo_.stochastic_steady_state = sss(oo_.dr, irfburninlength, options_.order);

for ivariable = 1:M_.endo_nbr //
      oo_.sss.(deblank(M_.endo_names(ivariable,:))) = oo_.stochastic_steady_state(ivariable,:);
end //

for ivariable = 1:M_.endo_nbr //
      oo_.ss.(deblank(M_.endo_names(ivariable,:))) = oo_.dr.ys(ivariable,:);
end //

for ishock = 1:M_.exo_nbr //
    variables_irfsss  =  irfsss(oo_.dr, M_.Sigma_e(:,ishock), irflength, irfburninlength, options_.order)';

    for ivariable = 1:M_.endo_nbr //
        oo_.irfs.(  strcat(deblank(M_.endo_names(ivariable,:)), '_', deblank(M_.exo_names(ishock,:)))  ) = variables_irfsss(:,ivariable);
    end //

end //



//-------------------------------------------------------------------------
//               Plot results
//-------------------------------------------------------------------------

varstoplot = {'U'            'pi'        'i'                     'v'          'J'            'p'};
varnames   = {'Unemployment' 'Inflation' 'Nominal interest rate' 'Vacancies'  'Match value'  'Relative price'};

level_vars = {};
irf_len = 20;
x_tick = 5;
subplot_x = 2;
subplot_y = 3;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);
