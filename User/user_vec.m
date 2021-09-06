% user_vec.m
%
% Define useful vectors for computing the risk-adjusted DSGE model
% solution. The elements of the eta vectors are filled in by user_mod.m.
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 26-Mar-2018

%%

eta_r_w = zeros(1, V.nmodel);
eta_c   = zeros(1, V.nmodel);
eta_r   = zeros(1, V.nmodel);
eta_r_h = zeros(1, V.nmodel);
eta_u   = zeros(1, V.nmodel);
eta_z   = zeros(1, V.nmodel);

H = [ para( P.vphi_u )^2 ; para( P.vphi_n )^2 ; 1 ];

%% END