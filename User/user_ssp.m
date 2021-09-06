% user_ssp.m
%
% Define pre-specified/implied parameters, nonstochastic steady-state
% ratios, and other entities. The file uses the parameter structure P and
% vector of parameter values para.
%
% Output
%
% ssp  vector of parameter values
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 26-Mar-2018

%%
% This section provides values for pre-specified parameters (ie,
% parameters that are not estimated).

ssp( P.alpha ) = 0.36;

ssp( P.delta ) = 0.025;

ssp( P.scale ) = 0.2;

%%
% This section transforms prior parameters (where necessary) before they
% appear in the canonical form of the DSGE model.

ssp( P.gamma ) = 1 + para( P.gamma_Q )/100;

ssp( P.beta ) = 1/(para( P.r_Q )/100 + 1);

%%
% This section defines steady-state variables and parameters used in the
% canonical form of the DSGE model.

ssp( P.ky ) = ssp( P.alpha )*ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi ))/(1 - (1 - ssp( P.delta ))...
    *ssp( P.beta )*ssp( P.gamma )^(-1/para( P.psi )));

ssp( P.iy ) = (1 - (1 - ssp( P.delta ))/ssp( P.gamma ))*ssp( P.ky );

ssp( P.wy ) = 1 - ssp( P.alpha );

ssp( P.dy ) = 1 - ssp( P.wy ) - ssp( P.iy );

ssp( P.cy ) = 1 - ssp( P.iy );

ssp( P.jy ) = ssp( P.cy )/(ssp( P.beta )...
    *ssp( P.gamma )^(1 - 1/para( P.psi )))/(1/(ssp( P.beta )...
    *ssp( P.gamma )^(1-1/para( P.psi ))) - 1);

ssp( P.hy ) = ssp( P.jy ) - ssp( P.wy ) - ssp( P.ky ) - ssp( P.dy );

%%
% This section defines other entities for streamlining the model equations.

ssp( P.zeta ) = (1 - para( P.phi ))/(1 - 1/para( P.psi ));

%% END