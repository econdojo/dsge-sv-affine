function logpdf = mvt_pdf(x, mu, Sigma, v)

% mvt_pdf.m
%
% Evaluate log multivariate Student-t density; uses common mean vector
% and scaling matrix as opposed to MATLAB built-in function mvnpdf.
%
% Input
%
% x      parameter row vectors (N x dim)
% mu     mean vector (1 x dim)
% Sigma  scaling matrix (dim x dim); covariance = v/(v-2)*Sigma
% v      degrees of freedom
%
% Output
%
% logpdf  log multivariate Student-t density (N x 1)
%
% Written by Dave Rapach and Fei Tan, Saint Louis University
%
% Updated: 17-Feb-2018

%%
% Initialization

[N, dim] = size(x);

Mu = repmat(mu, N, 1);

R = cholmod(Sigma);

Sigma = R'*R;

%%
% Choose density type

if isinf(v)  % multivariate normal

    const = -dim/2*log(2*pi) - log(abs(det(Sigma)))/2;

    logpdf = const - sum(((x-Mu)/Sigma).*(x-Mu), 2)/2;

else  % multivariate Student-t

    const = v/2*log(v) + log(gamma((dim + v)/2)) - dim/2*log(pi) ...
        - log(gamma(v/2)) - log(abs(det(Sigma)))/2;

    logpdf = const - (dim + v)/2 ...
        *log(abs(v + sum(((x-Mu)/Sigma).*(x-Mu), 2)));
end

%% END