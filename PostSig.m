function s2 = PostSig(s2, P, V, Y, pdof)

% PostSig.M
%
% Draw the posterior stochastic variance series using Gibbs sampling with
% a Metropolis-Hastings step and tailoring of the proposal density. We
% make draws for sigma2(0), sigma2(1), ..., sigma2(T). The stochastic
% variances are drawn sequentially based on the following conditional
% posterior kernel density:
%
% p[ sigma2(t) | sigma2(\t), theta, Y ]
%
% proportional to
%
% p[y(t) | sigma2(t-1), sigma2(t), Y(1:t-1), theta]*...
% p[y(t+1) | sigma2(t), sigma2(t+1), Y(1:t), theta]*...
% ...
% p[y(T) | sigma2(T - 1), sigma2(T), Y(1:T), theta]*...
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
% s2    [1 x (T + 1)] vector of (current) stochastic variances
% P     structure of model parameters
% V     structure of model variables
% Y     (T x number of observables) data matrix
% pdof  student-t degrees of freedom for proposal density
%
% Output
%
% s2  [1 x (T + 1)] vector of (new) stochastic variances
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% This section solves the model without risk adjustment. The solution is
% computed using gensys.m; user_mod.m defines the vectors/matrices serving
% as inputs for gensys.m, while user_vec.m and user_ssp.m define
% convenient vectors and steady-state ratios, respectively. Note that
% SSR.Sigma_e is a column vector containing the steady-state variances for
% the complete set of shocks (including the volatility shock); SSR.Sigma_u
% is a column vector containing the measurment error variances. The
% vectors/matrices for the state-space representation are stored in the
% structure SSR.

para = P.para( : , 1 );  % need for user_vec/ssp/mod/mod_affine

G0  = zeros(V.nmodel);
G1  = zeros(V.nmodel);
CC  = zeros(V.nmodel, 1);
Psi = zeros(V.nmodel, V.nshock);
Pi  = zeros(V.nmodel, V.nfore);

SSR.Sigma_e = zeros(V.nshock, 1);
SSR.Sigma_u = zeros(V.ndata, 1);

SSR.Z = zeros(V.ndata, V.nmodel);
SSR.D = zeros(V.ndata, 1);

user_vec
user_ssp

j = 0;  % need for user_mod

user_mod

[G_old, C_old, M_old, ~, ~, ~, ~, ~] = gensys(G0, G1, CC, Psi, Pi);

%%
% This section re-solves the model to obtain the risk-adjusted solution.
% We iterate until C, G, and M of the state-space representation converge;
% user_mod_affine.m specifies additional components needed to compute the
% affine form of the DSGE model allowing for risk adjustment. The affine
% form provides the inputs to gensys.m to compute the affine state-space
% representation of the risk-adjusted model solution.


for iter = 1:1000

    j = 0;  % need for user_mod_affine

    user_mod_affine

    [SSR.G, SSR.C, SSR.M, ~, ~, ~, ~, eu] = gensys(G0, G1, CC, Psi, Pi);
    
    if eu( 1 ) ~= 1 || eu( 2 ) ~= 1 || ...
            norm([ SSR.G SSR.C SSR.M ]-[ G_old C_old M_old ], Inf) < 1e-5

        break

    else

        G_old = SSR.G;
        C_old = SSR.C;
        M_old = SSR.M;

    end

end

%%
% This section computes matrices for the state transition and measurement
% equations for the reduced state vector. It also computes the ergodic
% moments for the reduced state vector.

rho_s = SSR.G( end , end );

Sigma_e = SSR.Sigma_e;

C1 = SSR.C( 1:end-1 );
C2 = SSR.C( end );

G12 = SSR.G( 1:end-1 , end );
M12 = SSR.M( 1:end-1 , end );
G11 = SSR.G( 1:end-1 , 1:end-1 );
M11 = SSR.M( 1:end-1 , 1:end-1 );

D  = SSR.D;
Z1 = SSR.Z( : , 1:end-1 );

Sigma_u = SSR.Sigma_u;

g0 = (eye(V.nmodel-1) - G11)...
    \(C1 - C2*M12 + (G12 + (1 - rho_s)*M12)*Sigma_e( end-1 )/100);

P_g0 = dlyap_sym(G11, M11*diag(Sigma_e( 1:end-1 ))*M11', 1);

%%
% This section simulates draws of sigma2(0),...,sigma2(T-1). The sigma2(0)
% draw is from the ergodic distribution. We impose the non-negativity
% constraint for the stochastic variance by truncating from below.

T = length(s2) - 1;

s2( 1 ) = Sigma_e( end-1 )/100 + sqrt(Sigma_e( end )/(1 - rho_s^2))*randn;

if s2( 1 ) < 0

    s2( 1 ) = 1e-4;

end

g_tm1 = g0;
P_tm1 = P_g0;

for t = 1:T-1

    s2_t = s2( t+1 );

    s2_extra = [ s2( t ) s2( t+2:end ) ];

    % csminwel.m is used to find the mode and Hessian of the log posterior
    % kernel density and tailor the mean vector and covariance matrix for
    % the proposal density. The negative of the log posterior kernel
    % density is computed by PostKer2.m.

    [~, mu_s2_t, ~, Hess_t, ~, ~, ~] = csminwel('PostKer2', ...
        s2_t, 1e-4, [], 1e-5, 100, 1, s2_extra, g_tm1, P_tm1, SSR, ...
        P, V, Y( t:end , : ), 0);

    % The Metropolis-Hastings step makes a draw for sigma2(t) from the
    % tailored student-t distribution. The acceptance probability is
    % computed by MoveProb.m.

    s2_t = mvt_rnd_mex(mu_s2_t, Hess_t, pdof, 1);

    [~, ~, alpha] = MoveProb('PostKer2', s2( t+1 ), s2_t, mu_s2_t, ...
        Hess_t, pdof, [], [], s2_extra, g_tm1, P_tm1, SSR, ...
        P, V, Y( t:end , : ), 0);

    if rand < alpha

        s2( t+1 ) = s2_t;

    end

    if s2( t+1 ) < 0

        s2( t+1 ) = 1e-4;

    end

    % The filtered mean vector and covariance matrix for the reduced state
    % vector for period t are computed by running another recursion of the
    % conditional Kalman filter. These moments are needed to subsequently
    % draw sigma2(t+1).

    s2_t = s2( t+1 );

    s2_tm1 = s2( t );

    C_t = C1 - C2*M12 + M12*s2_t + (G12 - rho_s*M12)*s2_tm1;

    Sigma_v_t = H*100*s2_tm1;

    g_predict = C_t + G11*g_tm1;
    P_predict = G11*P_tm1*G11' + M11*diag(Sigma_v_t)*M11';

    y_predict = D + Z1*g_predict;

    F = Z1*P_predict*Z1' + diag(Sigma_u);
    K = (P_predict*Z1')/F;

    g_tm1 = g_predict + K*(Y( t , : )' - y_predict);
    P_tm1 = P_predict - K*Z1*P_predict;

end

%%
% The final section draws sigma2(T). In this case, we ignore the
% p[sigma2(t + 1) | ...] component in the posterior kernel density.

s2_t = s2( end );

s2_extra = s2( end-1 );

[~, mu_s2_t, ~, Hess_t, ~, ~, ~] = csminwel('PostKer2', ...
    s2_t, 1e-2, [], 1e-5, 100, 1, s2_extra, g_tm1, P_tm1, SSR, ...
    P, V , Y( end , : ), 1);

s2_t = mvt_rnd_mex(mu_s2_t, Hess_t, pdof, 1);

[~, ~, alpha] = MoveProb('PostKer2', s2( end ), s2_t, mu_s2_t, Hess_t, ...
    pdof, [], [], s2_extra, g_tm1, P_tm1, SSR, P, V, Y( end , : ), 1);

if rand < alpha

    s2( end ) = s2_t;

end

if s2( end ) < 0

    s2( end ) = 1e-4;

end

%% END