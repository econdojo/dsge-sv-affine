function [P, V] = ParVar

% ParVar.m
%
% Construct parameter/variable structures
%
% Input
%
% none
%
% Output
%
% P  parameter structure
% V  variable structure
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

%%
% This section creates the parameter structure.          

user_parvar  % defines para/ssp/mvar/mshock/mfore/dvar

P.name  = para( : , 1 );                        % parameter names
P.type  = para( : , 5 );                        % parameter types
P.para  = cell2mat(para( : , [ 2:4 , 6:8 ] ));  % parameter values
P.npara = length(P.name);                       % # of model parameters
P.nssp  = length(ssp);                          % # of steady-state pars

% Parameter indices

for k = 1:P.npara

    eval([ 'P.' P.name{ k } ' = k;' ]);

end

% Steady-state indices

for k = 1:P.nssp

    eval([ 'P.' ssp{ k } ' = k;' ]);

end

% Keep parameters away from bounds (change if necessary)

P.realsmall = 1e-5;

%%
% This section creates the variable structure.          

V.mvar   = mvar;            % model variables
V.mshock = mshock;          % model structural shocks
V.nmodel = length(mvar);    % number of model variables
V.nshock = length(mshock);  % number of structural shocks
V.nfore  = length(mfore);   % number of forecast errors
V.ndata  = length(dvar);    % number of data variables

% Model variable indices

for k = 1:V.nmodel

    eval([ 'V.model.' mvar{ k } ' = k;' ]);

end

% Structural shock indices

for k = 1:V.nshock

    eval([ 'V.shock.' mshock{ k } ' = k;' ]);

end

% Forecast error indices

for k = 1:V.nfore

    eval([ 'V.fore.' mfore{ k } ' = k;' ]);

end

% Data variables indices

for k = 1:V.ndata

    eval([ 'V.data.' dvar{ k } ' = k;' ]);

end

%% END