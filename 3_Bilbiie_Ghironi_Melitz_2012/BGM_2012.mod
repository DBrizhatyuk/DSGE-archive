//-------------------------------------------------------------------------
// Bilbiie Ghironi Melitz (2012) Endogenous Entry, Product Variety, and Business Cycles
//-------------------------------------------------------------------------

var p   // varietry effect
    mu  // markup
    w   // wage
    d   // profit
    D   // aggregate profit
    c   // consumption
    N   // number of operating firms
    N_e // number of entrants
    I_E // investment in entry
    L   // labor
    L_E // labor in product creation
    L_C // labor in production
    v   // firm value
    R_E // ex ante return on product creation
    Z   // TFP
    y   // firm output
    Y   // aggegate output

    log_c log_p log_d log_v log_L log_w log_N log_Z log_N_e log_mu log_y log_Y log_D log_R_E log_I_E log_L_E log_L_C; //

varexo e_Z;    // TFP shock

parameters betta   //
           delta  //
           theta  //
           f_e    //
           chi    //
           phi    //
           N_bar sigm //
           rho_z; //




//-------------------------------------------------------------------------
//               Parameters
//-------------------------------------------------------------------------

betta = 0.99;            // discount factor
delta = 0.025;          // exit rate
theta = 3.8;            // elastitity of substitution b\w varieties
f_e = 1;                // entry cost
chi = 0.924271;         // weight of disutility of labor period utility
phi = 4;                // inverse elasticity of labor supply

rho_z = 0.979;          // TFP shock persistence

// translog case if 1, CES otherwise
@#define translog = 1

// translog preference paramaters
sigm = 0.35323;
N_bar = 1;




//-------------------------------------------------------------------------
//               Model equilibrium conditions (table 2)
//-------------------------------------------------------------------------

model;

// (1) Pricing
p = mu * w/Z;

// (2-3) Markup & Variety effect (update the steady-state file when choosing one or the other)
@#if translog == 0
   mu = 1 + 1/(sigm*N);
   p = exp( -0.5*(N_bar-N)/(sigm*N_bar*N) );
@#else
   mu = theta/(theta-1);
   p = N^(1/(theta-1));
@#endif

// (4) Profits
d = (1-1/mu) * c/N;

// (5) Free entry
v = w * f_e/Z;

// (6) Number of firms
N = (1-delta)*(N(-1) + N_e(-1));

// (7) Intratemporal optimality
chi*L^(1/phi) = w/c;

// (8) Euler equation (shares)
v = betta * (1-delta) * (c/c(+1)) * (v(+1) + d(+1));

// (9) Aggregate accounting
c + N_e*v = w*L + N*d;

// (10) Productitity
ln(Z) = rho_z*ln(Z(-1)) + e_Z;

//-------------------------------------------------------------------------

// (11) Firm output
y = c/(N*p);

// (12) Aggregate output
Y = c + N_e*v;

// (13) Aggregate profit
D = d*N;

// (14) Return to investment
R_E = (v(+1) + d(+1))/v;

// (15) Investment
I_E = N_e*v;

// (16) Labor in product creation
L_E = N_e*f_e/Z;

// (17) Labor in production
L_C = L - L_E;

// log varibles
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
var e_Z = 0.01^2;
end;

stoch_simul(irf = 100, order=1, nomoments, nograph, nofunctions);


varstoplot = {'log_c'        'log_Y' 'log_N_e'  'log_N'                'log_d'        'log_v'       'log_R_E'               'log_I_E'     'log_w' 'log_L'         'log_L_E'                   'log_L_C'             'log_p'           'log_y'        'log_D'             'log_Z'};
varnames   = {'Consumption'  'GDP'   'Entry'    'Number of products'   'Firm profit'  'Firm value'  'Return to investment'  'Investment'  'Wage'  'Total labor'  'Labor in product creation'  'Labor in production' 'Relative price'  'Firm output'  'Aggregate profit'  'TFP'};

level_vars = {};
irf_len = 100;
x_tick = 20;
subplot_x = 4;
subplot_y = 4;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);




//-------------------------------------------------------------------------
//               Stochastic simulations
//-------------------------------------------------------------------------

shocks;
var e_Z = 0.0072^2;
end;

stoch_simul(irf=0, order=1, periods=100000, nomoments, nograph, nofunctions);
