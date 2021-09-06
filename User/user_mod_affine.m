% user_mod_affine.m
%
% Fill in elements of CC and G0 to allow for computation of the risk-
% adjusted model solution. The file uses the parameter vector para and
% structure SSR.
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 08-Mar-2018

%%
% [1] Euler equation - aggregate wealth

j = j + 1;

rvec = (eta_r_w - 1/para( P.psi )*eta_c ...
    - ssp( P.beta )/(1 - ssp( P.beta ))*eta_u ...
    + (1 - 1/para( P.psi ))*eta_z)*M_old;

G0( j , V.model.sigma2 ) = ssp( P.zeta )*rvec*diag([ H ; 0 ])*rvec'/2;

CC( j ) = -ssp( P.zeta )...
    *rvec*diag([ zeros(length(H), 1) ; 1 ].*SSR.Sigma_e)*rvec'/200;

%%
% [2] Euler equation - equity

j = j + 1;

rvec = ((ssp( P.zeta ) - 1)*eta_r_w ...
    - 1/para( P.psi )*ssp( P.zeta )*eta_c ...
    + eta_r - ssp( P.zeta )*ssp( P.beta )/(1 - ssp( P.beta ))*eta_u ...
    + (1 - para( P.phi ))*eta_z)*M_old;

G0( j , V.model.sigma2 ) = rvec*diag([ H ; 0 ])*rvec'/2;

CC( j ) = -rvec*diag([ zeros(length(H), 1) ; 1 ].*SSR.Sigma_e)*rvec'/200;

%%
% [3] Risk-free return (non-detrended)

j = j + 1;

rvec = ((ssp( P.zeta ) - 1)*eta_r_w ...
    - 1/para( P.psi )*ssp( P.zeta )*eta_c ...
    - ssp( P.zeta )*ssp( P.beta )/(1 - ssp( P.beta ))*eta_u ...
    - para( P.phi )*eta_z)*M_old;

G0( j , V.model.sigma2 ) = rvec*diag([ H ; 0 ])*rvec'/2;

CC( j ) = -rvec*diag([ zeros(length(H), 1) ; 1 ].*SSR.Sigma_e)*rvec'/200;

%% END