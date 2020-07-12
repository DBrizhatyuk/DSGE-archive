function m = plot_irfs_dynare(varstoplot, varnames, level_vars, irf_len, x_tick, subplot_x, subplot_y)

%-------------------------------------------------------------------------
% Plots and saves a custom figure with Dynare IRFs
%
% INPUTS:
%
% varstoplot -- variables to plot
% varnames   -- titles of IRFs
% level_vars -- plot IRFs in levels instead of log deviations
% irf_len    -- number of IRF periods
% x_tick     -- x-axis labels every *x_tick* periods
% subplot_x  -- number of rows of subplots
% subplot_y  -- number of columns of subplots
%-------------------------------------------------------------------------

global M_ oo_
% M_         -- Dynare model setup
% oo_        -- Dynare simulation results

vars_number = length(varstoplot); % total number of variables to plot
steady_st   = zeros(1,irf_len);
x_axis = [0:irf_len];


for k=1:M_.exo_nbr      % Iterate over shocks
%figure('units','normalized','outerposition',[0 0 1 1])

  for i=1:vars_number   % Iterate over variables
    subplot(subplot_y,subplot_x,i)

    if any(strcmp(level_vars,varstoplot(i)))
      hold on
      plot(x_axis, [0 eval( char(strcat( 'oo_.irfs.', varstoplot(i), '_' , deblank(M_.exo_names(k,:)) ))) ], 'color', [0, 0.4470, 0.7410], 'LineWidth', 3)
      plot(steady_st, 'k--', 'LineWidth', 0.5)
      hold off
      ylabel('\Delta level', 'fontsize',12)

    else
      hold on
      plot(x_axis, [0 eval( char(strcat( 'oo_.irfs.', varstoplot(i), '_' , deblank(M_.exo_names(k,:)) ))) ]*100, 'color', [0, 0.4470, 0.7410], 'LineWidth', 3)
      plot(steady_st, 'k--', 'LineWidth', 0.5)
      hold off
      ylabel('% deviation', 'fontsize',12)

    end

	  box on
    grid on
    set(gca,'Xtick',0:x_tick:irf_len, 'fontsize',12)

	  xlim([0 irf_len])
	  title(varnames(i), 'fontsize',15)
    %xlabel('quarters')

  end


set(gcf, 'PaperSize', [4*subplot_x 3*subplot_y]);
set(gcf, 'PaperPosition', [0 0 4*subplot_x 3*subplot_y]);

fig_name = [deblank(M_.exo_names(k,:)) '.pdf'];
print('-dpdf','-r100', strcat(fig_name));
close;


end
end
