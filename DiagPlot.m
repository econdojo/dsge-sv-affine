function DiagPlot(P, chain)

% DiagPlot.m
%
% Plot convergence diagnostics
%
% Input
%
% P      parameter structure
% chain  MCMC posterior draws
%
% Output
%
% none
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% This section takes care of preliminaries for plotting.

M       = size(chain, 1);                   % sample size after burn-in's
step    = round([ M/50 (M/10):(M/10):M ]);  % step length
rm      = zeros(11, P.npara);               % recursive means
plotpos = find(P.para( : , 6 ));            % determine plot position
nplot   = length(plotpos);                  % number of plot parameters

for k = 1:11

    rm( k , : ) = mean(chain( 1:step( k ) , : ));

end

if nplot <= 4

    dim1 = 1;

elseif nplot >= 5 && nplot <= 8

    dim1 = 2;

elseif nplot >= 9 && nplot <= 12

    dim1 = 3;

elseif nplot >= 13

    dim1 = 4;

end

dim2 = ceil(nplot/dim1);

set(0, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultAxesFontSize', 12)

%%
% This section plots the recursive means.

figure(1)

for k = 1:nplot

    subplot(dim1, dim2, k);

    hold on
    grid on
    box on;

    plot(step, rm( : , plotpos( k ) ), 'r', 'linewidth', 1.5);

    title(strtrim(P.name{ plotpos( k ) }), 'fontsize', 12, ...
        'fontname', 'times');

end

saveas(gcf, 'User/Fig_ParMean.fig');

%%
% This section plots the autocorrelation functions.

figure(2)

for k = 1:nplot

    subplot(dim1, dim2, k);

    autocorr(chain( : , plotpos( k )), 500);

    title(strtrim(P.name{ plotpos( k ) }), 'fontsize', 12, ...
        'fontname', 'times');

end

saveas(gcf, 'User/Fig_ACF.fig');

%%
% This section plots the traceplots.

figure(3)

for k = 1:nplot

    subplot(dim1, dim2, k);

    hold on
    box on

    plot(1:M, chain( : , plotpos( k ) ), 'linewidth', 1.5);

    title(strtrim(P.name{ plotpos( k ) }), 'fontsize', 12, ...
        'fontname', 'times');

end

saveas(gcf, 'User/Fig_Trace.fig');

%% END