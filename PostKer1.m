function postker = PostKer1(theta, block, s2, P, V, Y)

% PostKer1.m
%
% Evaluate the (negative of the log) posterior kernel density taking the
% following form:
%
% p(theta | {sigma2}, Y)
%
% proportional to
%
% p(Y | {sigma2}, theta)*p({sigma2} | theta)*p(theta)
%
% To compute p(Y | {sigma2}, theta), we use a conditional Kalman filter
% based on the following system for the reduced state vector:
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
% theta  block of model parameters
% block  block of parameter indices
% s2     [1 x (T+1)] vector of stochastic variances
% P      structure of model parameters
% V      structure of model variables
% Y      (T x number of observables) data matrix
%
% Output
%
% postker  (negative of) log posterior kernel density
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% This section constructs the state-space representation of the model
% solution. We first check that the model parameters fall within the
% admissable parameter space.

if max(theta' < P.para( block , 2 )+P.realsmall) || ...
    max(theta' > P.para( block , 3 )-P.realsmall) 

    postker = 1e10;

    return

else

    % The model is first solved without risk adjustment. The solution is
    % computed using gensys.m; user_mod.m defines the vectors/matrices
    % serving as inputs for gensys.m, while user_vec.m and user_ssp.m
    % define conveninent vectors and steady-state parameters/ratios,
    % respectively. Note that SSR.Sigma_e is a column vector containing
    % the steady-state variances for the complete set of shocks (including
    % the volatility shock); SSR.Sigma_u is a column vector containing the
    % measurment error variances. The vectors/matrices for the state-space
    % representation are in the structure SSR. We restrict the parameter
    % space to unique saddle-path solutions.

    P.para( block , 1 ) = theta';

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

    [G_old, C_old, M_old, ~, ~, ~, ~, eu] = gensys(G0, G1, CC, Psi, Pi);

    if eu( 1 ) ~= 1 || eu( 2 ) ~= 1

        postker = 1e10;

        return

    end

end

% The models is then re-solved to obtain the risk-adjusted solution. We
% iterate until C, G, and M of the state-space representation converge;
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
% After checking that the parameter vector is consistent with a unique
% saddle-path solution, this section computes the negative of the log
% posterior kernel density.

if eu( 1 ) ~= 1 || eu( 2 ) ~= 1

    postker = 1e10;

    return

else

    % Vectors/matrices are computed for the state transition and
    % mesurement equations corresponding to the reduced state vector. The
    % intercept vector and covariance matrix for the shocks in the state
    % transition equation are time varying. Note that SSR.Sigma_e is
    % originally a column vector containing the steady-state variances for
    % the complete set of shocks (including the volatility shock).
    % SSR.Sigma_e is then redefined as a matrix, where each of its T
    % columns contains the variances for the non-volatility shocks for a
    % given period. Similarly, SSR.C is redefined from a column vector to
    % a matrix, where each of its T columns contains the intercept vector
    % for a given period for the reduced state vector.

    Sigma_e = SSR.Sigma_e;

    rho_s = SSR.G( end , end );

    T = length(s2) - 1;

    C1 = SSR.C( 1:end-1 );
    C2 = SSR.C( end );

    G12 = SSR.G( 1:end-1 , end );
    M12 = SSR.M( 1:end-1 , end );

    SSR.C = repmat(C1-C2*M12, 1, T) + M12*s2( 2:end ) ...
        + (G12 - rho_s*M12)*s2( 1:end-1 );

    SSR.G = SSR.G( 1:end-1 , 1:end-1 );
    SSR.M = SSR.M( 1:end-1 , 1:end-1 );

    SSR.Sigma_e = (repmat(H, 1, T))...
        .*(repmat(100*s2( 1:end-1 ), length(H), 1));

    SSR.Z = SSR.Z( : , 1:end-1 );

    g0 = (eye(V.nmodel-1) - SSR.G)...
        \(C1 - C2*M12 + (G12 + (1 - rho_s)*M12)*Sigma_e( end-1 )/100);

    P_g0 = dlyap_sym(SSR.G, SSR.M*diag(Sigma_e( 1:end-1 ))*SSR.M', 1);

    % The redefined form of SSR serves as an input for KalmanFilter.m,
    % which uses the conditional Kalman filter to compute the log
    % likelihood for Y. The log likelihhood for {sigma2} is then computed
    % using mvt_pdf.m. Finally, prior.m is used to compute the log
    % marginal prior densities to complete the log posterior kernel
    % density.
    
    [~, ~, loglik_Y] = KalmanFilter_mex(Y, SSR, g0, P_g0);

    eps_s = s2( 2:end ) - (1 - rho_s)*Sigma_e( end-1 )/100 ...
        - rho_s*s2( 1:end-1 );

    loglik_s2 = mvt_pdf(eps_s', 0, Sigma_e( end ), Inf);

    logprior = 0;

    for k = 1:P.npara

        logprior = logprior ...
            + Prior(P.para( k , 1 ), P.para( k , 4 ), P.para( k , 5 ), ...
            P.type{ k });
    end

    postker = -(sum(loglik_Y) + sum(loglik_s2) + logprior);

end

%% END