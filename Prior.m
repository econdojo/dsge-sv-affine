function logprior = Prior(x, para1, para2, type)

% Prior.m
%
% Evaluate log prior density
%
% Input
%
% x      parameter value
% para1  prior parameter 1
% para2  prior parameter 2
% type   prior type
%
% Output
%
% logprior  log prior density
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%% 

switch type

    case 'G'  % Gamma distribution (para1 = mean, para2 = std dev)

        a = para1^2/para2^2;
        b = para2^2/para1;

        logprior = log(gampdf(x, a, b));

    case 'N'  % Normal distribution (para1 = mean, para2 = std dev)

        logprior = log(normpdf(x, para1, para2));

    case 'B'  % Beta distribution (para1 = mean, para2 = std dev)

        a = -para1*(para2^2 + para1^2 - para1)/para2^2;
        b = (para1 - 1)*(para2^2 + para1^2 - para1)/para2^2;

        logprior = log(betapdf(x, a, b));

    case 'V'  % Inverse-Gamma distribution (para1 = nu, para2 = s2)

        a = para1/2;
        b = a*para2;

        logprior = log((b^a/gamma(a))*x^(-(a + 1))*exp(-b/x));

    case 'I'  % Inverse-Gamma type-I distribution (para1 = nu, para2 = s)

        a = para1/2;
        b = a*para2^2;

        logprior = log((2*b^a/gamma(a))*x^(-(2*a + 1))*exp(-b/x^2));

    case 'U'  % Uniform dist (para1 = lower bound, para2 = upper bound)

        logprior = log(unifpdf(x, para1, para2));

    otherwise

        error('Prior type does not exist!')

end

%% END