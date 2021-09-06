function [s_smooth, P_smooth] = KalmanSmoother(SSR, s_filter, P_filter)

% KalmanSmoother.m
%
% Compute smoothed densities. The Kalman filter is based on the
% following system:
%
% s(t) = C + G*s(t-1) + M*epsilon(t)
% E_t[epsilon(t)*epsilon(t)'] = Sigma_e
% y(t) = D + Z*s(t) + u(t)
% E[u(t)*u(t)'] = Sigma_u
%
% Input
%
% SSR       state-space representation (structure)
% s_filter  state filtering mean from Kalman filter
% P_filter  state filtering covariance from Kalman filter
%
% Output
%
% s_smooth(:,t)  E[s(t)|Y(:,1:T)] for t = 0,...,T
% P_smooth(:,t)  Cov[s(t)|Y(:,1:T)] for t = 0,...,T
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%% 
% This section handles preliminaries for Kalman smoothing.

[dims, T] = size(s_filter);

s_smooth = [ zeros(dims, T-1) s_filter( : , end ) ];
P_smooth = cat(3, zeros(dims, dims, T-1), P_filter( : , : , end ));


%%
% This section performs the Kalman smoothing recursions.

for t = (T-1):-1:1

    ps = SSR.C( : , t ) + SSR.G*s_filter( : , t );

    Omega_ps = SSR.G*P_filter( : , : , t )*SSR.G' ...
        + SSR.M*diag(SSR.Sigma_e( : , t ))*SSR.M';
    
    gain = (P_filter( : , : , t )*SSR.G')*pinv(Omega_ps);

    s_smooth( : , t ) = s_filter( : , t ) ...
        + gain*(s_filter( : , t+1 ) - ps);

    P_smooth( : , : , t ) = P_filter( : , : , t ) ...
        + gain*(P_filter( : , : , t+1 ) - Omega_ps)*gain';
end

%% END