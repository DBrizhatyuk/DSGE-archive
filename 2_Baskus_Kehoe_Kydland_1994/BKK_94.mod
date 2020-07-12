//--------------------------------------------------------------------------------------------------------------
// Bakus Kehoe Kydland (1994): Dynamics of the Trade Balance and the Terms of Trade: The J-Curve?
//--------------------------------------------------------------------------------------------------------------




var f_H f_F // 2 output
    y_H y_F // 4 aggregate basket
    c_H c_F // 6 consumption
    n_H n_F // 8 labor
    l_H l_F // 10 leisure
    w_H w_F // 12 wage
    k_H k_F // 14 capital
    I_H I_F // 16 investment
    q_H q_F // 18 Tobin's q
    R_k_H R_k_F // 20 Return on capital

    nx_H nx_F     // 22 net exports
    MU_c_H MU_c_F // 24 marginal utility of consumption
    sdf_H sdf_F   // 26 stochastic discount factor

    qa_H qa_F // 28 relative price of the home good at home and foreign markets
    qb_H qb_F // 30 relative price of the foreign good at home and foreign markets
    a_H a_F   // 32 home good at home and foreign markets
    b_H b_F   // 34 foreign good at home and foreign markets
    tot       // 35 terms of trade

    A_H A_F   // 37 TFP
    log_y_H log_y_F log_c_H log_c_F log_I_H log_I_F log_k_H log_k_F log_n_H log_n_F log_A_H log_A_F;

varexo u_H u_F; // TFP shocks

parameters betta gama muu      //
           theta delta psi_k   //
           omega sigm          //
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
psi_k = 0;       // capital adjustment cost

sigm   = 1.5;     // elasticity of substitution between home and foreign goods
omega  = 0.85;    // home bias in consumption basket

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
f_@{X} = exp(A_@{X}) * k_@{X}(-1)^(theta)*n_@{X}^(1-theta);

// (3-4) Labor supply`
((1-muu)/muu)*(c_@{X}/l_@{X}) = w_@{X};

// (5-6) Capital supply`
q_@{X} = sdf_@{X} * ( (1-delta)*q_@{X}(+1) + R_k_@{X}(+1) );

// (7-8) Tobin's q
q_@{X}*( 1 - (psi_k/2)*(I_@{X}/I_@{X}(-1)-1)^2 - psi_k*(I_@{X}/I_@{X}(-1))*(I_@{X}/I_@{X}(-1)-1)  )  + sdf_@{X}*q_@{X}(+1)*psi_k*(I_@{X}(+1)/I_@{X}-1)*(I_@{X}(+1)/I_@{X})^2 = 1;

// (9-10) Capital accumulation
k_@{X} = (1-delta)*k_@{X}(-1) + (1 - (psi_k/2)*(I_@{X}/I_@{X}(-1)-1)^2)*I_@{X};

// (11-12) Marginal utility of consumption
MU_c_@{X} = muu * (c_@{X}^B1) * (l_@{X}^B2);

// (13-14) Stochastic discount factor
sdf_@{X} = betta*MU_c_@{X}(+1)/MU_c_@{X};

// (15-16) Labor - Leisure
n_@{X} = 1 - l_@{X};

// (17-18) Aggregate good resourse constraints
y_@{X} = I_@{X} + c_@{X};

@#endfor

// (19-20) Labor demand`
w_H = qa_H * (1-theta)*exp(A_H)  * (k_H(-1)/n_H)^theta;
w_F = qb_F * (1-theta)*exp(A_F)  * (k_F(-1)/n_F)^theta;

// (21-22) Capital demand
R_k_H = qa_H * theta*exp(A_H) * (k_H(-1)/n_H)^(theta-1);
R_k_F = qb_F * theta*exp(A_F) * (k_F(-1)/n_F)^(theta-1);

// (23-24) Home and foreign goods resourse constraints
f_H = a_H + a_F;
f_F = b_F + b_H;

// (25-26) Domestic good demands
a_H = omega*(qa_H)^(-sigm) * y_H;
b_F = omega*(qb_F)^(-sigm) * y_F;

// (27-28) Imported good demands
b_H = (1-omega)*(qb_H)^(-sigm) * y_H;
a_F = (1-omega)*(qa_F)^(-sigm) * y_F;

// (29-30) Aggregate price index
omega*(qa_H)^(1-sigm) + (1-omega)*(qb_H)^(1-sigm) = 1;
omega*(qb_F)^(1-sigm) + (1-omega)*(qa_F)^(1-sigm) = 1;

// (31-32) Net exports to GDP
nx_H = ( a_F - b_H*(qb_H/qa_H) )/y_H;
nx_F = ( b_H - a_F*(qa_F/qb_F) )/y_F;

// (33-34) Perfect risk-sharing
MU_c_H/MU_c_F = qa_F/qa_H;
MU_c_H/MU_c_F = qb_F/qb_H;

// (35) Term of trade: foreign/domestic
tot = qb_H/qa_H;

// (36-37) Productivity process
A_H = (1-rho1-rho2) + rho1*A_H(-1) + rho2*A_F(-1) + u_H;
A_F = (1-rho1-rho2) + rho1*A_F(-1) + rho2*A_H(-1) + u_F;

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


varstoplot = {'log_y_H' 'log_c_H' 'log_I_H' 'nx_H'         'log_A_H'  'tot'                    'log_y_F'    'log_c_F'    'log_I_F'    'nx_F'            'log_A_F'};
varnames   = {'Y home'  'C home'  'I home'  'NX/GDP home'  'TFP home' 'Terms of trade'         'Y foreign'  'C foreign'  'I foreign'  'NX/GDP foreign'  'TFP foreign'};

level_vars = {'nx_H', 'nx_F', 'tot'};
irf_len = 60;
x_tick = 10;
subplot_x = 6;
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

stoch_simul(irf=0, order=1, periods=100000, nomoments, nograph, nofunctions) log_y_H log_y_F log_c_H log_c_F log_I_H log_I_F log_k_H log_k_F log_n_H log_n_F nx_H tot;


[ytrend_H, ycyclical_H] = sample_hp_filter( [log_y_H log_c_H log_I_H log_k_H log_n_H, nx_H, tot], 1600 );
[ytrend_F, ycyclical_F] = sample_hp_filter( [log_y_F log_c_F log_I_F log_k_F log_n_F], 1600 );

standard_devs_H         = std([ycyclical_H])*100;

fprintf('Standard deviations:\n');
fprintf('%20s \t %3.2f\n','sigma(y)',standard_devs_H(1));
fprintf('%20s \t %3.2f\n','sigma(c)',standard_devs_H(2));
fprintf('%20s \t %3.2f\n','sigma(I)',standard_devs_H(3));
fprintf('%20s \t %3.2f\n','sigma(k)',standard_devs_H(4));
fprintf('%20s \t %3.2f\n','sigma(n)',standard_devs_H(5));
fprintf('%20s \t %3.2f\n','sigma(nx/y)',standard_devs_H(6));
fprintf('%20s \t %3.2f\n\n','sigma(tot)',standard_devs_H(7));

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
fprintf('%20s \t %3.2f\n','corr(y_H, nx_H)',corr(ycyclical_H(:,1), ycyclical_H(:,6)));
fprintf('%20s \t %3.2f\n','corr(y_H, tot)',corr(ycyclical_H(:,1), ycyclical_H(:,7)));
fprintf('%20s \t %3.2f\n','corr(nx_H, tot)',corr(nx_H, ycyclical_H(:,7)));


//-------------------------------------------------------------------------

J_curve = [];
k=1;

for i=1:11
  J_curve(1,k) = corr( ycyclical_H(6:(end-5),7), ycyclical_H(i:(end-11+i),6) );
  k = k+1;
end

lag = [-5:5];

figure
hold on
plot(lag, J_curve, 'LineWidth', 3)
plot(lag, zeros(1,11), 'k--', 'LineWidth', 0.5)
plot(zeros(1,10), [-1:0.25:1, 1], 'k--', 'LineWidth', 0.5)
box on
grid on

set(gca,'Xtick',-5:1:5, 'fontsize',12)
title('J curve: corr(tot_{t}, nx_{t+j})', 'fontsize',15)
xlabel('lag, j')
ylim([-0.61 0.61])



set(gcf, 'PaperSize', [8 6]);
set(gcf, 'PaperPosition', [0 0 8 6]);

print('-dpdf','-r100', 'J_curve.pdf');
close;
