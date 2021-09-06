function e_hat = DistSmoother(Y, SSR, s0, P0)

% DistSmoother.m
%
% Compute smoothed shocks; see Koopman (1993), 'Disturbance smoother for
% state space models', Biometrika. The smoothed shocks are computed using
% output from the Kalman filter, which is based on the following system:
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
% e_hat(:,t)  E[e(t)|Y(:,1:T)] for t = 1,...,T
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%% 
% This section handles preliminaries for Kalman filtering.

T = size(Y, 1);  % number of periods

[dim_y, dim_s] = size(SSR.Z);  % number of observables/states

dim_e = size(SSR.Sigma_e, 1);  % number of shocks

u_hat = zeros(dim_y, T);  % prediction error

F = zeros(dim_y, dim_y, T);  % prediction covariance
K = zeros(dim_s, dim_y, T);  % Kalman gain

e_hat = zeros(dim_e, T);  % smoothed shocks

s_filter = s0;  % filtered state mean
P_filter = P0;  % filtered state covariance

%%
% This section performs the Kalman recursions.

for t = 1:T

    % State prediction uses the state-transition equation to compute the
    % mean vector and covariance matrix for the one-step-ahead prediction
    % of the state vector.

    s_predict = SSR.C( : , t ) + SSR.G*s_filter;

    P_predict = SSR.G*P_filter*SSR.G' ...
        + SSR.M*diag(SSR.Sigma_e( : , t ))*SSR.M';
    
    % Observable prediction uses the measurment equation to compute the
    % mean vector and covariance matrix for the one-step-ahead prediction
    % of the vector of observable variables.

    y_hat = SSR.D + SSR.Z*s_predict;

    u_hat( : , t ) = Y( t , : )' - y_hat;

    F( : , : , t ) = SSR.Z*P_predict*SSR.Z' + diag(SSR.Sigma_u);

    % Updating provides the filtered states via the Kalman gain.

    K( : , : , t ) = (P_predict*SSR.Z')/F( : , : , t );

    s_filter = s_predict + K( : , : , t)*u_hat( : , t );
    P_filter = P_predict - K( : , : , t)*SSR.Z*P_predict;

end

%%
% This section proceeds in reverse using output from the Kalman filter
% to compute the smoothed disturbances.

r = zeros(dim_s, 1);

for t = T:-1:1

    r = (SSR.Z'/F( : , : , t ))*u_hat( : , t ) ...
        + (SSR.G - SSR.G*K( : , : , t )*SSR.Z)'*r;

    e_hat( : , t ) = diag(SSR.Sigma_e( : , t ))*SSR.M'*r;

end

%% END