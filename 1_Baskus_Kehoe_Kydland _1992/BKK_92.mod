//--------------------------------------------------------------------------------------------------------------
// Bakus Kehoe Kydland (1992): International Real Business Cycles
//--------------------------------------------------------------------------------------------------------------




var y_H y_F // 2 output
    c_H c_F // 4 consumption
    n_H n_F // 6 labor
    l_H l_F // 8 leisure
    w_H w_F // 10 wage
    k_H k_F // 12 capital
    I_H I_F // 14 investment
    q_H q_F // 16 Tobin's q
    R_k_H R_k_F // 18 Return on capital

    nx_H nx_F     // 20 net exports
    MU_c_H MU_c_F // 22 marginal utility of consumption
    sdf_H sdf_F   // 24 stochastic discount factor

    A_H A_F   // 26 TFP

    log_y_H log_y_F log_c_H log_c_F log_I_H log_I_F log_k_H log_k_F log_n_H log_n_F log_A_H log_A_F;

varexo u_H u_F; // TFP shocks

parameters betta gama muu      //
           theta delta psi_k   //
           rho1 rho2           //
           B1 B2;              //




//-------------------------------------------------------------------------
//               Parameters
//-------------------------------------------------------------------------

betta = 0.99;    // discount factor
gama = -1;       // (1-gamma) = elasticity of intertemporal sub.; risk-aversion
muu   = 0.34;    // share of consumption in utility (Cobb-Douglas utility)

theta = 0.36;    // elasticity of output wrt capital
delta = 0.025;   // depreciation rate
psi_k = 0.02;    // capital adjustment cost

rho1 = 0.906;    // persistence on the productivity shock
rho2 = 0.088;    // cross-country producticity shock spillover

B1 = muu*(gama-1) + (muu-1);
B2 = (1-muu)*(gama-1) + (1-muu);




//-------------------------------------------------------------------------
//               Model equilibrium conditions
//-------------------------------------------------------------------------

model;

@#for X in [ "H", "F" ]

// (1-2) Production function
y_@{X} = exp(A_@{X}) * k_@{X}(-1)^(theta)*n_@{X}^(1-theta);

// (3-4) Labor supply`
((1-muu)/muu)*(c_@{X}/l_@{X}) = w_@{X};

// (5-6) Labor demand`
w_@{X} = (1-theta)*exp(A_@{X})  * (k_@{X}(-1)/n_@{X})^theta;

// (7-8) Capital supply`
q_@{X} = sdf_@{X} * ( (1-delta)*q_@{X}(+1) + R_k_@{X}(+1) );

// (9-10) Tobin's q
q_@{X}*( 1 - (psi_k/2)*(I_@{X}/I_@{X}(-1)-1)^2 - psi_k*(I_@{X}/I_@{X}(-1))*(I_@{X}/I_@{X}(-1)-1)  )  + sdf_@{X}*q_@{X}(+1)*psi_k*(I_@{X}(+1)/I_@{X}-1)*(I_@{X}(+1)/I_@{X})^2 = 1;

// (11-12) Capital accumulation
k_@{X} = (1-delta)*k_@{X}(-1) + (1 - (psi_k/2)*(I_@{X}/I_@{X}(-1)-1)^2)*I_@{X};

// (13-14) Capital demand
R_k_@{X} = theta*exp(A_@{X}) * (k_@{X}(-1)/n_@{X})^(theta-1);

// (15-16) Marginal utility of consumption
MU_c_@{X} = muu * (c_@{X}^B1) * (l_@{X}^B2);

// (17-18) Stochastic discount factor
sdf_@{X} = betta*MU_c_@{X}(+1)/MU_c_@{X};

// (19-20) Labor - Leisure
n_@{X} = 1 - l_@{X};

// (21-22) Net exports to GDP
nx_@{X} = (y_@{X} - c_@{X} - I_@{X})/y_@{X};

@#endfor

// (23-26) Productivity process
A_H = (1-rho1-rho2) + rho1*A_H(-1) + rho2*A_F(-1) + u_H;
A_F = (1-rho1-rho2) + rho1*A_F(-1) + rho2*A_H(-1) + u_F;


// (27) World resource constraint
y_H + y_F = (c_H + c_F) + (I_H + I_F);

// (28) Perfect risk-sharing
MU_c_H =  MU_c_F;

// (29) - (38) Log variables
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
var u_H = 0.00852^2;
var u_F = 0.00852^2;
corr u_H, u_F = 0;
end;

stoch_simul(irf = 60, order=1, nomoments, nograph, nofunctions);


varstoplot = {'log_y_H' 'log_c_H' 'log_I_H' 'nx_H'         'log_A_H'          'log_y_F'    'log_c_F'    'log_I_F'    'nx_F'            'log_A_F'};
varnames   = {'Y home'  'C home'  'I home'  'NX/GDP home'  'TFP home'         'Y foreign'  'C foreign'  'I foreign'  'NX/GDP foreign'  'TFP foreign'};

level_vars = {'nx_H', 'nx_F'};
irf_len = 60;
x_tick = 10;
subplot_x = 5;
subplot_y = 2;

plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y);




//-------------------------------------------------------------------------
//               Stochastic simulations
//-------------------------------------------------------------------------

shocks;
var u_H = 0.00852^2;
var u_F = 0.00852^2;
corr u_H, u_F = 0.258;
end;

stoch_simul(irf=0, order=1, periods=100000, nomoments, nograph, nofunctions) log_y_H log_c_H log_c_F log_y_F;


[ytrend_H, ycyclical_H] = sample_hp_filter( [log_y_H log_c_H log_I_H log_k_H log_n_H], 1600 );
[ytrend_F, ycyclical_F] = sample_hp_filter( [log_y_F log_c_F log_I_F log_k_F log_n_F], 1600 );

standard_devs_H          = std([ycyclical_H, nx_H])*100;

fprintf('Standard deviations:\n');
fprintf('%20s \t %3.2f\n','sigma(y)',standard_devs_H(1));
fprintf('%20s \t %3.2f\n','sigma(c)',standard_devs_H(2));
fprintf('%20s \t %3.2f\n','sigma(I)',standard_devs_H(3));
fprintf('%20s \t %3.2f\n','sigma(k)',standard_devs_H(4));
fprintf('%20s \t %3.2f\n','sigma(n)',standard_devs_H(5));
fprintf('%20s \t %3.2f\n\n','sigma(nx/y)',standard_devs_H(6));

relative_standard_devs_H = standard_devs_H./standard_devs_H(1);

fprintf('Relative standard deviations:\n');
fprintf('%20s \t %3.2f\n','sigma(y)',relative_standard_devs_H(1));
fprintf('%20s \t %3.2f\n','sigma(c)',relative_standard_devs_H(2));
fprintf('%20s \t %3.2f\n','sigma(I)',relative_standard_devs_H(3));
fprintf('%20s \t %3.2f\n','sigma(k)',relative_standard_devs_H(4));
fprintf('%20s \t %3.2f\n\n','sigma(n)',relative_standard_devs_H(5));

country_corr = corr( ycyclical_H, ycyclical_F );

fprintf('Correlations:\n');
fprintf('%20s \t %3.2f\n','corr(y_H, y_F)',country_corr(1,1));
fprintf('%20s \t %3.2f\n','corr(c_H, c_F)',country_corr(2,2));
fprintf('%20s \t %3.2f\n','corr(y_H, nx_H)',corr(ycyclical_H(:,1), nx_H));
