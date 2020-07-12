//--------------------------------------------------------------------------------------------------------------
// Gertler and Karadi (2011) A model of unconventional monetary policy
//--------------------------------------------------------------------------------------------------------------




var C     // Consumption
    W     // Real wage
    L     // Labor
    R_K   // Retutn on capital
    U     // Capital utilization
    delta_K // Capital depreciation
    mu    //
    Omega //
    nu    //
    N     //
    S     //
    K     // Effective capital
    B     //
    Y     // Output
    R     // Real intrest rate
    P     // Ralative price = 1/markup
    pi    // Inflation
    I     // Investment
    I_n   //
    q     // Tobin's q
    spread leverage //
    lambda sdf sdf_adjusted  //
    log_xi log_Y log_C log_I log_K log_L log_q log_N //
    xi     // Capital quality
    i;     // Nominal interest rate

varexo e_xi; // Capital quality shock

parameters betta delta_K_ss delta_B gama psi alphaa theta epsilon omega pi_ss i_ss chi psi_k psi_p h rho_i kappa_pi kappa_y Y_ss I_ss c1 c2 A //
           rho_xi;




//-------------------------------------------------------------------------
//               Parameters
//-------------------------------------------------------------------------

betta    = 0.99;    // Discount factor
h        = 0.815;   // Habit persistence
gama     = 1.0001;  // Inverse of elasticity of intertemporal substitution
psi      = (0.276)^(-1);  // Elasticity of labor supply

delta_K_ss  = 0.025;   // Capital depreciation rate
psi_k     = 1.728;    // Capital adjustment scale
alphaa    = 0.33;    // Capital share

theta    = 0.381;   // Fraction of capital that can be diverted
delta_B  = 0.972;   // Survival rate of the bankers
omega    = 0.002;   // Proportional transfer to the entering bankers

epsilon  = 4.16;    // Elasticity of substitution b\w retailers varieties
psi_p    = 120;     // Price adjustment costs (Rottemberg)
pi_ss    = 1;       // Steady state inflation
Y_ss     = 1;       // Steady state output
i_ss     = ( betta )^(-1) * pi_ss - 1;  // Steady state home nominal interest rate
rho_i    = 0.8;     // Smoothing parameter of the Taylor rule
kappa_y  = 0.50/4;  // Output gap coef. of the Taylor rule
kappa_pi = 1.5;     // Inflation coef. of the Taylor rule

c2 = 7.2;

rho_xi  = 0.6; // Capital quality shock persistence

// Initial values for the steady-state solution
chi   = 1;
A     = 1;
I_ss  = 1;
c1    = 1;

@#define GHH = 0




///-------------------------------------------------------------------------
//               Model equilibrium conditions
//-------------------------------------------------------------------------

model;

// (1) Final good market clearing
Y = C + I + (psi_k/2)*( (I_n+I_ss)/(I_n+I_ss) - 1 )^2 * (I_n+I_ss);

//  (2) Final output
Y = A * (U*K(-1)*xi)^(alphaa) * L^(1-alphaa);

// (3) Labor supply
@#if GHH == 0
    W = chi * L^(1/psi) / lambda;
@#else
    W = chi * L^(1/psi);
@#endif

// (4) Labor demand
W = P*(1-alphaa)*(Y/L);

// (5-6) Euler equation (deposits)
sdf * R(+1) = 1;
R = (1+i(-1))/pi;

// (7-8) Capital accumulation
K = xi*K(-1) + I_n;
I_n = I - delta_K*xi*K;

// (9) Capital utilization
P*alphaa*Y/U = (c1 + c2*(U-1)) * xi*K(-1);

// (10) Capital depreciation
delta_K = delta_K_ss + c1*(U-1) + (c2/2)*(U-1)^2;

// (11) Capital demand
R_K = (P*alphaa*Y/(xi*K(-1)) + (1-delta_K)*q)*xi / q(-1);

// (12) Tobin's q
q = 1 + (psi_k/2)*( (I_n+I_ss)/(I_n(-1)+I_ss) - 1 )^2 + psi_k*( (I_n+I_ss)/(I_n(-1)+I_ss) )*( (I_n+I_ss)/(I_n(-1)+I_ss) - 1 ) - psi_k*sdf*( (I_n(+1)+I_ss)/(I_n+I_ss) )^2 * ( (I_n(+1)+I_ss)/(I_n+I_ss) - 1 );

// (13) Time-varying markup
1/P = epsilon/( (epsilon-1) + psi_p*(pi/pi_ss)*(pi/pi_ss-1) - psi_p*sdf*(pi(+1)/pi_ss-1)*(pi(+1)/pi_ss)*(Y(+1)/Y) );

// (14) Taylor rule
i = (1-rho_i)*( i_ss + kappa_pi*(pi-pi_ss) + kappa_y*( epsilon/(epsilon-1) - 1/P ) ) + rho_i*i(-1);

// (15) Spread
theta * mu = sdf * Omega(+1) * (R_K(+1) - R(+1));

// (16) Ex ante MV of net worth
Omega = (1-delta_B) + delta_B*nu;

// (17) Ex post MV of net worth
nu = (sdf*Omega(+1)*R(+1))/(1-mu);

// (18) Net worth
N = delta_B * ( (R_K-R)*S(-1)*q(-1) + R*N(-1) ) + omega*S(-1)*q(-1);

// (19) Binding constraint
nu*N = theta*S*q;

// (20) Banks resourse constraint
B + N = S*q;

// (21)
K = S;

// (22) Leverage
leverage = S*q/N;

//-------------------------------------------------------------------------

// (23-24)
sdf = betta * lambda(+1)/lambda;
sdf_adjusted =  betta * lambda(+1)/lambda * Omega(+1);

// (25)
@#if GHH == 1
    lambda = ( C - h*C(-1) - chi*L^(1+1/psi)/(1+1/psi) )^(-gama) - h*( C(+1) - h*C - chi*L(+1)^(1+1/psi)/(1+1/psi) )^(-gama);
@#else
    lambda = (C - h*C(-1))^(-gama) - betta*h*(C(+1) - h*C)^(-gama);
@#endif

// (26)
spread = R_K(+1) - R;

//-------------------------------------------------------------------------

// (27) Capital quality shock
log(xi) = rho_xi*log(xi(-1)) - e_xi;

//-------------------------------------------------------------------------

log_xi = log(xi);
log_Y = log(Y);
log_C = log(C);
log_I = log(I);
log_K = log(K);
log_L = log(L);
log_q = log(q);
log_N = log(N);

end;




//-------------------------------------------------------------------------
//               Steady state
//-------------------------------------------------------------------------

steady;
resid(1);

// Save steady-state values for convenience
load('steady_st_values.mat');
NumberOfEndogenousVariables = M_.orig_endo_nbr;
for ii = 1:NumberOfEndogenousVariables
  oo_.steady_state_vars(1).( strcat( deblank(M_.endo_names(ii,:)) )) = steady_st_values(ii);
end


//-------------------------------------------------------------------------
//               Steady state
//-------------------------------------------------------------------------

steady (maxit = 10000, solve_algo = 2);
check;




//-------------------------------------------------------------------------
//               Impulse responses
//-------------------------------------------------------------------------

shocks;
var e_xi = 0.05^2;
end;

stoch_simul(irf = 40, order=1, nomoments, nograph, nofunctions);


varstoplot = {'log_xi' 'R' 'spread'      'log_Y' 'log_C' 'log_I' 'log_K' 'log_L' 'log_q' 'log_N'  'pi' 'i'};
varnames   = {'\xi'    'R' 'E(R_{K})-R'  'Y'     'C'     'I'     'K'     'L'     'Q'     'N'     '\pi' 'i'};

level_vars = {'R', 'spread', 'pi', 'i'};
irf_len = 40;
x_tick = 10;
subplot_x = 3;
subplot_y = 4;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);
