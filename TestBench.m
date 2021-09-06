% TestBench.m
%
% Macro-finance DSGE model analysis & generate MEX functions
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

clear
clc
addpath('User', 'Toolbox', 'MacMex', 'WinMex');

%%
% Solve model

[P, V] = ParVar;

Y = importdata('Data.txt');

para = P.para( : , 1 );

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

j = 0;

user_mod

[G_old, C_old, M_old, ~, ~, ~, ~, eu] = gensys(G0, G1, CC, Psi, Pi);

for iter = 1:1000

    j = 0;

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

fprintf('Existence  =  %d,  Uniqueness  =  %d\n', eu( 1 ), eu( 2 ));

%%
% Unconditional Kalman filter/smoother

T = length(Y);

C = SSR.C;

Sigma_e = SSR.Sigma_e;

SSR.C = repmat(C, 1, T);

SSR.Sigma_e = repmat(Sigma_e, 1, T);

s0 = (eye(V.nmodel) - SSR.G)\C;
P0 = dlyap_sym(SSR.G, SSR.M*diag(Sigma_e)*SSR.M', 1);

[sf, Pf, ~] = KalmanFilter_mex(Y, SSR, s0, P0);

[ss, Ps] = KalmanSmoother_mex(SSR, sf, Pf);

es = DistSmoother_mex(Y, SSR, s0, P0);

%%
% Conditional Kalman filter

s2 = mean(chain_s2);

rho_s = SSR.G( end , end );

C1 = C( 1:end-1 );
C2 = C( end );

G12 = SSR.G( 1:end-1 , end );
M12 = SSR.M( 1:end-1 , end );

SSR.C = repmat(C1-C2*M12, 1, T) + M12*s2( 2:end ) ...
    + (G12 - rho_s*M12)*s2( 1:end-1 );

SSR.G = SSR.G( 1:end-1 , 1:end-1 );
SSR.M = SSR.M( 1:end-1 , 1:end-1 );

SSR.Sigma_e = (repmat(H, 1, T)).*(repmat(100*s2( 1:end-1 ), length(H), 1));

SSR.Z = SSR.Z( : , 1:end-1 );

g0 = (eye(V.nmodel-1) - SSR.G)...
    \(C1 - C2*M12 + (G12 + (1 - rho_s)*M12)*Sigma_e( end-1 )/100);

P_g0 = dlyap_sym(SSR.G, SSR.M*diag(Sigma_e( 1:end-1 ))*SSR.M', 1);

[g_filter, P_filter, loglik] = KalmanFilter_mex(Y, SSR, g0, P_g0);

%% END