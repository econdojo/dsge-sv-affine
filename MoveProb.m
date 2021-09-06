function [pk_last, pk_next, alpha] = MoveProb(fun, x, y, mu, Sigma, ...
    pdof, pk_last, pk_next, varargin)

% MoveProb.m
%
% Compute Metropolis-Hastings probability of move; see Chib & Greenberg
% (1995), 'Understanding the Metropolis-Hastings algorithm', American
% Statistician. The proposal density is a multivariate t-distribution.
%
% Input
%
% fun       string naming posterior kernel function
% x         current parameter vector (1 x dim)
% y         candidate parameter vector (1 x dim)
% mu        proposal mean vector (1 x dim)
% Sigma     proposal scaling matrix (dim x dim)
% pdof      proposal t degrees of freedom
% pk_last   log current posterior kernel
% pk_next   log candidate posterior kernel
% varargin  additional inputs required by fun.m
%
% Output
%
% pk_last  log current posterior kernel
% pk_next  log candidate posterior kernel
% alpha    M-H probability of move
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% The log posterior kernel for the current parameter is computed (if
% necessary) using fun.m. Note that fun.m returns the negative of
% the log posterior kernel.

if isempty(pk_last)

    pk_last = -feval(fun, x, varargin{ : });

end

%%
% The log posterior kernel for the candidate parameter is computed (if
% necessary) using fun.m. Note that fun.m returns the negative of
% the log posterior kernel.

if isempty(pk_next)

    pk_next = -feval(fun, y, varargin{ : });

end

%%
% The probability of a move is computed using the current and candidate
% log posterior kernel and proposal density values.

r = exp(pk_next - pk_last + mvt_pdf_mex(x, mu, Sigma, pdof) ...
    - mvt_pdf_mex(y, mu, Sigma, pdof));

alpha = min([ r 1 ]);

%% END