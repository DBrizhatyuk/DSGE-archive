//-----------------------------------------------------------------------------------
// Jesús Fernández-Villaverde et al. (2011): Risk Matters: The Real Effects of Volatility Shocks
//-----------------------------------------------------------------------------------




var C // consumption
    D // net foreing debt
    Y // output
    L // labor
    I // investment
    K // capital
    q // Tobin's q
    NX // trade balance
    CA // current account
    lambda // MU of C
    R e_R sigma_R; // interest rate process

varexo u_R u_sigma_R; // interest rate level and volatility shocks

parameters beta alph gamma psi chi delta phi_K A //
           R_ss D_ss phi rho_R sigma_R_ss rho_sigma_R sigma_sigma_R; // interest rate process

// 1 for GHH case   0 for the separable power utility
@#define GHH = 0




//-------------------------------------------------------------------------
//               Parameters (Argentina calibration, M1)
//-------------------------------------------------------------------------

beta  = 0.980;       // Discount factor
R_ss  = 1/beta - 1;  // Steady state world interest rate
alph = 0.32;         // Capital share
gamma = 5;           // Inverse of intertemporal elasticity of substitution
psi   = 1000;        // Inverse elasticity of labor supply
delta = 0.014;       // Capital depreciation rate
phi   = 0.001;       // Sensitivity of the interest rate to debt variations
phi_K = 95;          // Scale for capital adjustment cost

D_ss  = 4;           // Steady state net foreign debt
A     = 1;           // TFP
chi   = 1;           // Share of Labor in period utility

// Interest rate process parameters
rho_R = 0.97;         // interest rate level shock persistence
sigma_R_ss = -5.71;   // interest rate level shock st. dev. is exp(sigma_R_ss)
rho_sigma_R = 0.94;   // interest rate uncertainty shock persistence
sigma_sigma_R = 0.46; // interest rate uncertainty shock st. dev.




//-------------------------------------------------------------------------
//               Model equilibrium conditions
//-------------------------------------------------------------------------

model;

// (1) Capital accumulation
K = (1-delta)*K(-1) + (1 - phi_K/2*(I/I(-1)-1)^2 )*I;

// (2) Resourse constraint
Y - C - I = D(-1) - D/(1+R) +  (phi/2)*(D-D_ss)^2;

// (3) Trade balance
NX = Y - C - I;

// (4) MU_C
@#if GHH == 1
  lambda = ( C - chi*L^(1+psi)/(1+psi) )^(-gamma);
@#else
  lambda = C^(-gamma);
@#endif

// (5) Consumption Euler
lambda/(1+R) = lambda*phi*(D-D_ss) +  beta*lambda(+1);

// (6) Capital mrkt equilibrium
q = beta * lambda(+1)/lambda * ( (1-delta)*q(+1) + alph*Y(+1)/K );

// (7) Labor mrkt equilibrium
@#if GHH == 1
  chi*L^(psi) = (1-alph) * Y/L;
@#else
  chi*L^(psi) / lambda = (1-alph) * Y/L;
@#endif

// (8) Tobin's q
1 = q * ( 1 - phi_K/2 * (I/I(-1)-1)^2 - phi_K*I/I(-1)*(I/I(-1)-1) ) + beta * (lambda(+1)/lambda) * q(+1)*phi_K*(I(+1)/I)^2 * (I(+1)/I - 1);

// (9) Production function
Y = A * L^(1-alph) * K(-1)^(alph);

// (10) Current account
CA = D - D(-1);

// (11-13) World real interest rate process
R       = R_ss + e_R;
e_R     = rho_R*e_R(-1) + exp(sigma_R)*u_R;
sigma_R = (1-rho_sigma_R)*sigma_R_ss + rho_sigma_R*sigma_R(-1) + sigma_sigma_R*u_sigma_R;

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
var u_R = 1;
var u_sigma_R = 1;
end;

stoch_simul(order=3,periods=0,irf=0,noprint,nograph,nomoments,nofunctions,nocorr,pruning);

irfburninlength   = 5000;
irflength         = 45;

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

varstoplot = {'C' 'I' 'Y' 'L' 'R' 'D'};
varnames   = {'C' 'I' 'Y' 'L' 'R' 'D'};

level_vars = {};
irf_len = 45;
x_tick = 15;
subplot_x = 6;
subplot_y = 1;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);
