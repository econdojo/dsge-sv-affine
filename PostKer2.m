function postker = PostKer2(s2_t, s2_extra, g_tm1, P_tm1, SSR, P, V, ...
    Y_r, terminal)

% PostKer2.m
%
% Evaluate the (negative of the log) posterior kernel density taking the
% following form:
%
% p[sigma2(t) | sigma2(\t), theta, Y]
%
% proportional to
%
% p[y(t) | sigma2(t-1), sigma2(t), Y(1:t-1), theta]*...
% p[y(t+1) | sigma2(t), sigma2(t+1), Y(1:t), theta]*...
% ...
% p[y(T) | sigma2(T-1), sigma2(T), Y(1:T), theta]*...
% p[sigma2(t) | sigma2(t-1), theta]*...
% p[sigma2(t+1) | sigma2(t), theta]
%
% To compute p[y(t) | ...] and related terms, we use a conditional
% Kalman filter based on the following system for the reduced state
% vector:
%
% g(t) = C(t) + G11*g(t-1) + M11*v(t), 
% C(t) = C1 - C2*M12 + M12*sigma2(t) + (G12 - rho*M12)*sigma2(t-1)
% E_t[v(t)*v(t)'] = H*sigma2(t-1)
% y(t) = D + Z1*g(t) + u(t)
% E[u(t)*u(t)'] = Sigma_u
%
% The reduced state vector excludes sigma2(t) from the complete state
% vector. The conditional Kalman filter derives from the following system
% for the complete state vector:
%
% s(t) = C + G*s(t-1) + M*epsilon(t)
% E_t[epsilon(t)*epsilon(t)'] = Sigma_e(t)
% y(t) = D + Z*s(t) + u(t)
% E[u(t)*u(t)'] = Sigma_u
%
% where s(t) = [ g(t)' sigma2(t) ]'.
%
% Stochastic volatility follows an AR(1) process:
%
% sigma2(t) = (1 - rho_s)*sigma2 + rho_s*sigma2(t-1) + epsilon_s(t)
%
% Input
%
% s2_t      sigma2(t)
% s2_extra  [ sigma2(t-1) sigma2(t+1) ... sigma2(T) ]
% g_tm1     t-1 filtered mean - reduced state vector
% P_tm1     t-1 filtered covariance matrix - reduced state vector
% SSR       structure of state-space representation elements
% P         structure of model parameters
% V         structure of model variables (needed for user_vec.m)
% Y_r       (T' x number of observables) data matrix for T' <= T
% terminal  (scalar) indicator: = 0 (1) for t < T (t = T)]
%
% Output
%
% postker  (negative of) log posterior kernel density
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% The following generates H (a column vector). Note that para is required
% for user_vec.m, which in turn provides H.

para = P.para( : , 1 );  % need for user_vec

user_vec

%%
% This section computes vectors/matrices for the state transition and
% measurement equations corresponding to the reduced state vector. The
% intercept vector and covariance matrix for the shocks in the state
% transition equation are time varying and correspond to period t. Note
% that SSR.Sigma_e is first a column vector containing the steady-state
% variances for the complete set of shocks (including the volatility
% shock); it is then redefined as a matrix with each column containing the
% vector of time-varying variances for a given period. SSR.Sigma_u is a
% column vector containing the measurment error variances.

if terminal == 1

    s2_r = [ s2_extra( 1 ) s2_t ];

elseif terminal == 0

    s2_r = [ s2_extra( 1 ) s2_t s2_extra( 2:end ) ];

end

Sigma_e = SSR.Sigma_e;

rho_s = SSR.G( end , end );

T_prime = length(s2_r) - 1;

C1 = SSR.C( 1:end-1 );
C2 = SSR.C( end );

G12 = SSR.G( 1:end-1 , end );
M12 = SSR.M( 1:end-1 , end );

SSR.C = repmat(C1-C2*M12, 1, T_prime) + M12*s2_r( 2:end ) ...
    + (G12 - rho_s*M12)*s2_r( 1:end-1 );

SSR.G = SSR.G( 1:end-1 , 1:end-1 );
SSR.M = SSR.M( 1:end-1 , 1:end-1 );

SSR.Sigma_e = repmat(H, 1, T_prime)...
    .*repmat(100*s2_r( 1:end-1 ), length(H), 1);

SSR.Z = SSR.Z( : , 1:end-1 );

%%
% This section computes the log likelihood for the observable variables
% for periods t through T. The log likelihood is computed using
% KalmanFilter based on the reduced state vector.

[~, ~, loglik] = KalmanFilter_mex(Y_r, SSR, g_tm1, P_tm1);

%%
% This section computes the conditional moments for sigma2(t) and
% sigma2(t + 1). It then computes the negative of the log posterior
% kernel density.

mu_s2_t = (1 - rho_s)*Sigma_e( end-1 )/100 + rho_s*s2_extra( 1 );

var_s2_t = Sigma_e( end );

if terminal == 1

    % For period T, we ignore p[sigma2(t + 1) | ...] in the posterior
    % kernel density.

    postker = -(sum(loglik) + mvt_pdf_mex(s2_t, mu_s2_t, var_s2_t, Inf));

elseif terminal == 0

    mu_s2_tp1 = (1 - rho_s)*Sigma_e( end-1 )/100 + rho_s*s2_t;

    var_s2_tp1 = Sigma_e( end );

    postker = -(sum(loglik) + mvt_pdf_mex(s2_t, mu_s2_t, var_s2_t, Inf) ...
        + mvt_pdf_mex(s2_extra( 2 ), mu_s2_tp1, var_s2_tp1, Inf));

end

%% END