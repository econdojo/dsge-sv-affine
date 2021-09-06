function [s_filter, P_filter, loglik] = KalmanFilter(Y, SSR, s0, P0)

% KalmanFilter.m
%
% Compute filtering densities & log likelihood. The Kalman filter is
% based on the following system:
%
% s(t) = C + G*s(t-1) + M*epsilon(t)
% E_t[epsilon(t)*epsilon(t)'] = Sigma_e
% y(t) = D + Z*s(t) + u(t)
% E[u(t)*u(t)'] = Sigma_u
%
% Input
%
% Y    (T x number of observables) data matrix
% SSR  state-space representation (structure)
% s0   initial filtering mean
% P0   initial filtering covariance
%
% Output
%
% s_filter(:,t)  E[s(t)|Y(:,1:t)] for t = 0,...,T
% P_filter(:,t)  Cov[s(t)|Y(:,1:t)] for t = 0,...,T
% loglik         log[p(y(t)|Y(:,1:t-1))] for t = 0,...,T
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%% 
% This section handles preliminaries for Kalman filtering.

T = size(Y, 1);  % number of periods

dim_s = length(s0);  % number of state variables

s_filter = [ s0 zeros(dim_s, T) ];              % filtering mean
P_filter = cat(3, P0, zeros(dim_s, dim_s, T));  % filtering covariance

loglik = zeros(T+1, 1);  % period log likelihood

%%
% This section performs the Kalman recursions.

for t = 1:T

    % State prediction uses the state-transition equation to compute the
    % mean vector and covariance matrix for the one-step-ahead prediction
    % of the state vector.

    s_predict = SSR.C( : , t ) + SSR.G*s_filter( : , t );

    P_predict = SSR.G*P_filter( : , : , t )*SSR.G' ...
        + SSR.M*diag(SSR.Sigma_e( : , t ))*SSR.M';
    
    % Observable prediction uses the measurment equation to compute the
    % mean vector and covariance matrix for the one-step-ahead prediction
    % of the vector of observable variables.

    y_hat = SSR.D + SSR.Z*s_predict;

    F = SSR.Z*P_predict*SSR.Z' + diag(SSR.Sigma_u);
    
    % The period log likelihood is based on the prediction-error
    % decomposition; mvt_pdf.m is used to compute the log likelihood.

    loglik( t+1 ) = mvt_pdf(Y( t , : ), y_hat', F, inf);

    % Updating provides the filtered states via the Kalman gain.

    K = (P_predict*SSR.Z')/F;

    s_filter( : , t+1 ) = s_predict + K*(Y( t , : )' - y_hat);

    P_filter( : , : , t+1 ) = P_predict - K*SSR.Z*P_predict;

end

%% END