function ineff = IneffFactor(chain, nlag)

% IneffFactor.m
%
% Compute the inefficiency factors of a MCMC sample based on the Parzen
% window: see Paolo Giordani, Michael Pitt, & Robert Kohn (2011), Bayesian
% inference for time series state space models, The Oxford Handbook of
% Bayesian Econometrics
%
% Input
%
% chain  univariate time series
% nlag   number of lags (<=200)
%
% Output
%
% ineff  inefficiency factor
%
% Note: This program is part of Dynare software
%
% Modified by Dave Rapach and Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%%
% Initialization

M = length(chain);

if nlag >= M

    error('Number of lags exceeds sample size minus one')

end

if mod(nlag, 2)

    nlag = nlag - 1;

end

%%
% Calculate autocorrelation function

acf = zeros(nlag+1, 1);

acf( 1 ) = 1;

m  = mean(chain);
sd = std(chain);

for k = 1:nlag

    acf( k+1 )=(chain( k+1:end ) - m)'*(chain( 1:end-k ) - m)./(M*sd^2);

end

%%
% Calculate Parzen weights

parzen = zeros(nlag+1, 1);

for k = 1:nlag/2+1

    parzen( k ) = 1 - 6*(k/nlag)^2 + 6*(k/nlag)^3;

end

for k = (nlag/2)+1:nlag+1

    parzen( k ) = 2*(1 - (k/nlag))^3;

end

parzen = parzen';

%%
% Calculate inefficiency factor

ineff = 1 + 2*sum(parzen( : ).*acf);

%% END