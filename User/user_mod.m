% user_mod.m
%
% Define matrices/vectors for the canonical form of the DSGE model and
% the measurement equation. The canonical form of the model is given by
%
% G0*s(t) = CC + G1*s(t-1) + Psi*epsilon(t) + Pi*eta(t),
% E[epsilon(t)*epsilon(t)'] = Sigma_e
%
% The measurement equation is given by
%
% y(t) = D + Z*s(t) + u(t)
% E[u(t)*u(t)'] = Sigma_u
%
% The file uses the parameter structure P, variable structure V, vector of
% parameter values para, and vector of values ssp.
%
% Output
%
% G0,G1,Psi,Pi,CC  vectors/matrices for GENSYS input
% Z,D              measurement equation vector/matrix
% Sigma_e          main diagonal of system covariance matrix
% Sigma_u          main diagonal of measurement covariance matrix
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 26-Mar-2018

%%
% This section fills in elements of CC, G0, G1, Psi, and Pi to define the
% DSGE model equations. It also fills in elements of the eta vectors
% used to compute the risk-adjusted solution.

% [1] Euler equation - aggregate wealth

j = j + 1;

G0( j , V.model.E_c ) = -1/para( P.psi );

G0( j , V.model.c ) = 1/para( P.psi );

G0( j , V.model.E_r_w ) = 1;

G0( j , V.model.upsilon ) = (1 - ssp( P.beta )*para( P.rho_u ))...
    /(1 - ssp( P.beta ));

G0( j , V.model.z ) = (1 - 1/para( P.psi ))*para( P.rho_z );

eta_r_w( V.model.r_w ) = 1;

eta_c( V.model.c ) = 1;

eta_u( V.model.upsilon ) = 1;

eta_z( V.model.z ) = 1;

% [2] Euler equation - equity

j = j + 1;

G0( j , V.model.E_c ) = -1/para( P.psi )*ssp( P.zeta );

G0( j , V.model.c ) = 1/para( P.psi )*ssp( P.zeta );

G0( j , V.model.E_r_w ) = ssp( P.zeta ) - 1;

G0( j , V.model.E_r ) = 1;

G0( j , V.model.upsilon ) = ssp( P.zeta )...
    *(1 - ssp( P.beta )*para( P.rho_u ))/(1 - ssp( P.beta ));

G0( j , V.model.z ) = (1 - para( P.phi ))*para( P.rho_z );

eta_r( V.model.r ) = 1;

% [3] Risk-free return (non-detrended)

j = j + 1;

G0( j , V.model.E_c ) = -1/para( P.psi )*ssp( P.zeta );

G0( j , V.model.c ) = 1/para( P.psi )*ssp( P.zeta );

G0( j , V.model.E_r_w ) = ssp( P.zeta ) - 1;

G0( j , V.model.r_f ) = 1;

G0( j , V.model.upsilon ) = ssp( P.zeta )...
    *(1 - ssp( P.beta )*para( P.rho_u ))/(1 - ssp( P.beta ));

G0( j , V.model.z ) = -para( P.phi )*para( P.rho_z );

% [4] Aggregate wealth return - definition

j = j + 1;

G0( j , V.model.r_w ) = 1;

G0( j , V.model.j ) = -1;

G1( j , V.model.j ) = -1/(ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi )));

G1( j , V.model.c ) = 1/(ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi ))) - 1;

% [5] Equity return (equals return to investment)

j = j + 1;

G0( j , V.model.r ) = 1;

G0( j , V.model.y ) = -(1 - (1 - ssp( P.delta ))*ssp( P.beta )...
    *ssp( P.gamma )^(-1/para( P.psi )));

G0( j , V.model.i ) = -ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi ))/para( P.xi )...
    *(1 - (1 - ssp( P.delta ))/ssp( P.gamma ));

G0( j , V.model.z ) = -(ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi ))/para( P.xi ) ...
    - (1 + 1/para( P.xi ))*(1 - ssp( P.delta ))*ssp( P.beta )...
    *ssp( P.gamma )^(-1/para( P.psi )));

G0( j , V.model.q ) = -(1 - ssp( P.delta ))*ssp( P.beta )...
    *ssp( P.gamma )^(-1/para( P.psi ));

G1( j , V.model.k ) = -(1 + ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi ))/para( P.xi ) ...
    - (1 + 1/para( P.xi ))*(1 - ssp( P.delta ))*ssp( P.beta )...
    *ssp( P.gamma )^(-1/para( P.psi )));

G1( j , V.model.q ) = -1;

% [6] Human wealth return - definition

j = j + 1;

G0( j , V.model.r_h ) = 1;

G0( j , V.model.h ) = -ssp( P.beta )*ssp( P.gamma )^(1 - 1/para( P.psi ));

G0( j , V.model.w ) = -(1 - ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi )));

G1( j , V.model.h ) = -1;

eta_r_h( V.model.r_h ) = 1;

% [7] Human wealth

j = j + 1;

G0( j , V.model.h ) = 1;

G0( j , V.model.j ) = -ssp( P.jy )...
    /(ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy ));

G0( j , V.model.w ) = ssp( P.wy )...
    /(ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy ));

G0( j , V.model.k ) = ssp( P.ky )...
    /(ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy ));

G0( j , V.model.q ) = ssp( P.ky )...
    /(ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy ));

G0( j , V.model.d ) = ssp( P.dy )...
    /(ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy ));

% [8] Marginal q

j = j + 1;

G0( j , V.model.q ) = 1;

G0( j , V.model.i ) = -1/para( P.xi );

G0( j , V.model.z ) = -1/para( P.xi );

G0( j , V.model.nu ) = 1;

G1( j , V.model.k ) = -1/para( P.xi );

% [9] Production function

j = j + 1;

G0( j , V.model.y ) = 1;

G0( j , V.model.z ) = ssp( P.alpha );

G1( j , V.model.k ) = ssp( P.alpha );

% [10] Capital accumulation

j = j + 1;

G0( j , V.model.k ) = 1;

G0( j , V.model.i ) = -(1 - (1 - ssp( P.delta ))/ssp( P.gamma ));

G0( j , V.model.nu ) = -(1 - (1 - ssp( P.delta ))/ssp( P.gamma ));

G0( j , V.model.z ) = (1 - ssp( P.delta ))/ssp( P.gamma );

G1( j , V.model.k ) = (1 - ssp( P.delta ))/ssp( P.gamma );

% [11] Wage

j = j + 1;

G0( j , V.model.w ) = 1;

G0( j , V.model.y ) = -1;

% [12] Dividend payment

j = j + 1;

G0( j , V.model.d ) = ssp( P.dy );

G0( j , V.model.y ) = -1;

G0( j , V.model.w ) = ssp( P.wy );

G0( j , V.model.i ) = ssp( P.iy );

% [13] Resource constraint

j = j + 1;

G0( j , V.model.y ) = 1;

G0( j , V.model.c ) = -(1 - ssp( P.iy ));

G0( j , V.model.i ) = -ssp( P.iy );

% [14] AR process - time-preference shock

j = j + 1;

G0( j , V.model.upsilon ) = 1;

G1( j , V.model.upsilon ) = para( P.rho_u );

Psi( j , V.shock.eps_u ) = 1;

% [15] AR process - investment shock

j = j + 1;

G0( j , V.model.nu ) = 1;

G1( j , V.model.nu ) = para( P.rho_n );

Psi( j , V.shock.eps_n ) = 1;

% [16] AR process - technology shock

j = j + 1;

G0( j , V.model.z ) = 1;

G1( j , V.model.z ) = para( P.rho_z );

Psi( j , V.shock.eps_z ) = 1;

% [17] AR process (w/intercept) - volatility shock

j = j + 1;

G0( j , V.model.sigma2 ) = 1;

G1( j , V.model.sigma2 ) = para( P.rho_s );

CC( j ) = (1 - para( P.rho_s ))*para( P.sigma )^2/100;

Psi( j , V.shock.eps_s ) = 1;

% [18] E_(t-1)[c(t)] error

j = j + 1;

G0( j , V.model.c ) = 1;

G1( j , V.model.E_c ) = 1;

Pi( j , V.fore.c ) = 1;

% [19] E_(t-1)[r_w(t)] error

j = j + 1;

G0( j , V.model.r_w ) = 1;

G1( j , V.model.E_r_w ) = 1;

Pi( j , V.fore.r_w ) = 1;

% [20] E_(t-1)[r(t)] error

j = j + 1;

G0( j , V.model.r ) = 1;

G1( j , V.model.E_r ) = 1;

Pi( j , V.fore.r ) = 1;

% [21] Lagged consumption

j = j + 1;

G0( j , V.model.d_c ) = 1;

G1( j , V.model.c ) = 1;

% [22] Lagged investment

j = j + 1;

G0( j , V.model.d_i ) = 1;

G1( j , V.model.i ) = 1;

% [23] Lagged wage

j = j + 1;

G0( j , V.model.d_w ) = 1;

G1( j , V.model.w ) = 1;

% [24] Lagged risk-free return

j = j + 1;

G0( j , V.model.d_r_f ) = 1;

G1( j , V.model.r_f ) = 1;

%%
% This section fills in elements of Z and D to define the measurment
% equations.

SSR.Z( V.data.EXR , V.model.r )     = para( P.vtheta );
SSR.Z( V.data.EXR , V.model.z )     = para( P.vtheta );
SSR.Z( V.data.EXR , V.model.d_r_f ) = -para( P.vtheta );

SSR.Z( V.data.RFR , V.model.d_r_f ) = 1;

SSR.Z( V.data.CGR , V.model.c )   = 1;
SSR.Z( V.data.CGR , V.model.z )   = 1;
SSR.Z( V.data.CGR , V.model.d_c ) = -1;

SSR.Z( V.data.IGR , V.model.i )   = 1;
SSR.Z( V.data.IGR , V.model.z )   = 1;
SSR.Z( V.data.IGR , V.model.d_i ) = -1;

SSR.Z( V.data.WGR , V.model.w )   = 1;
SSR.Z( V.data.WGR , V.model.z )   = 1;
SSR.Z( V.data.WGR , V.model.d_w ) = -1;

SSR.D( V.data.RFR ) = 100*log(ssp( P.gamma )^(1/para( P.psi ))...
    /ssp( P.beta ));

SSR.D( V.data.CGR ) = para( P.gamma_Q );

SSR.D( V.data.IGR ) = para( P.gamma_Q );

SSR.D( V.data.WGR ) = para( P.gamma_Q );

%%
% This section define the main diagonals of the covariance matrices for
% the DSGE model shocks and measurment errors.

SSR.Sigma_e( V.shock.eps_u ) = (para( P.vphi_u )*para( P.sigma ))^2;
SSR.Sigma_e( V.shock.eps_n ) = (para( P.vphi_n )*para( P.sigma ))^2;
SSR.Sigma_e( V.shock.eps_z ) = para( P.sigma )^2;
SSR.Sigma_e( V.shock.eps_s ) = (para( P.omega )/100)^2;

SSR.Sigma_u = (ssp( P.scale )*std(Y)').^2;

%% END