function x = mvt_rnd(mu, Sigma, v, N)

% mvt_rnd.m
%
% Generate multivariate Student-t distributed random vectors; uses
% uses common mean vector and scaling matrix as opposed to MATLAB
% built-in function mvnrnd.
%
% Input
%
% mu     mean vector (1 x dim)
% Sigma  scaling matrix (dim x dim); covariance = v/(v-2)*Sigma
% v      degrees of freedom
% N      number of realizations
%
% Output
%
% x  parameter row vectors (N x dim)
%
% Written by Dave Rapach and Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%%
% Initialization

dim = length(mu);

Mu = repmat(mu, N, 1);

R = cholmod(Sigma);

%%
% Choose distribution type

if isinf(v) % multivariate normal 

    x = Mu + randn(N, dim)*R;

else  % multivariate Student-t

    x = Mu + repmat(sqrt(v./chi2rnd(v, [ N 1 ])), 1, dim)...
        .*(randn(N, dim)*R);

end

%% END