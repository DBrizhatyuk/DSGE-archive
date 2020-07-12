function y = irfsss(dr, e1, long, drop, iorder)

% function y = irf(dr, e1, long, drop, replic, iorder)
% Computes impulse response functions from no-shock stochastic steady state
%
% INPUTS
%    dr:     structure of decisions rules for stochastic simulations
%    e1:     exogenous variables value in time 1 after one shock
%    long:   number of periods of simulation
%    drop:   truncation (in order 2)
%    iorder: first or second order approximation
%
% OUTPUTS
%    y:      impulse response matrix
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_ oo_ options_

temps = repmat(dr.ys,1,M_.maximum_lag);

ex1 = zeros(long+drop, M_.exo_nbr);

ex2 = ex1;
ex2(drop+1,:) = ex2(drop+1,:)+e1';


y1 = simult_(temps,dr,ex1,iorder);
y2 = simult_(temps,dr,ex2,iorder);

y = 0;

y = y + (   y2(:,M_.maximum_lag+drop+1:end) - y1(:,M_.maximum_lag+drop+1:end)  );

end
