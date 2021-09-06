% user_parvar.m
%
% Define DSGE model parameters and variables set by the user; ParVar.m
% loads this file as it creates parameter/variable structures for the DSGE
% model. The canonical form of the model is given by
%
% G0*s(t) = CC + G1*s(t-1) + Psi*epsilon(t) + Pi*eta(t)
%
% Output
%
% para    model parameters
% ssp     steady state/implied parameters
% mvar    model variables
% mshock  model structural shocks
% mfore   model forecast errors
% dvar    data variables
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 02-Apr-2018

%%
% This section specifies a matrix with information for the DSGE model
% parameters to be estimated. The columns contain the following
% information:
%
% Column 1 - parameter name
% Column 2 - initial value
% Column 3 - lower bound
% Column 4 - upper bound
% Column 5 - prior type as defined in Prior.m
% Column 6 - prior para(1) in Prior.m
% Column 7 - prior para(2) in Prior.m
% Column 8 - indicator to plot parameter

para = {
        'psi    ', 2.00, 0, 1e5, 'G', 2.00, 0.30, 1  % [1]
        'phi    ', 6.00, 0, 1e5, 'G', 6.00, 2.00, 1  % [2]
        'xi     ', 0.50, 0, 1e5, 'G', 1.50, 0.50, 1  % [3]
        'r_Q    ', 0.30, 0, 1e5, 'G', 0.30, 0.10, 1  % [4]
        'gamma_Q', 0.40, 0, 1e5, 'G', 0.40, 0.10, 1  % [5]
        'rho_u  ', 0.80, 0, 1,   'B', 0.80, 0.05, 1  % [6]
        'rho_n  ', 0.80, 0, 1,   'B', 0.80, 0.05, 1  % [7]
        'rho_z  ', 0.20, 0, 1,   'B', 0.20, 0.10, 1  % [8]
        'vphi_u ', 0.50, 0, 1e5, 'G', 1.00, 0.50, 1  % [9]
        'vphi_n ', 4.00, 0, 1e5, 'G', 1.00, 0.50, 1  % [10]
        'sigma  ', 0.70, 0, 1e5, 'I', 4.00, 0.40, 1  % [11]
        'rho_s  ', 0.80, 0, 1e5, 'B', 0.80, 0.05, 1  % [12]
        'omega  ', 0.10, 0, 1e5, 'I', 4.00, 0.05, 1  % [13]
        'vtheta ', 1.50, 0, 1e5, 'G', 1.50, 0.10, 1  % [14]
        };

%%
% This section defines names for pre-specified/implied parameters, non-
% stochastic steady-state ratios, and other entities in the DSGE model.

ssp = {
       'alpha'  % [1] pre-specified
       'delta'  % [2] pre-specified
       'scale'  % [3] pre-specified
       'gamma'  % [4] implied
       'beta'   % [5] implied
       'omega'  % [6] implied
       'ky'     % [7] K/Y nonstochastic steady-state ratio
       'iy'     % [8] I/Y nonstochastic steady-state ratio
       'wy'     % [9] W/Y nonstochastic steady-state ratio
       'dy'     % [10] D/Y nonstochastic steady-state ratio
       'cy'     % [11] C/Y nonstochastic steady-state ratio
       'jy'     % [12] J/Y nonstochastic steady-state ratio
       'hy'     % [13] H/Y nonstochastic steady-state ratio
       'zeta'   % [14] zeta = (1 - phi)/(1 - 1/psi)
       };

%%
% This section defines names for the model variables. The number of
% variables needs to equal the number of model equations defined in
% USER_MODEL.M.

mvar = {
        'y'       % [1] output
        'c'       % [2] consumption
        'i'       % [3] investment
        'k'       % [4] capital
        'w'       % [5] wage
        'j'       % [6] aggregate wealth
        'h'       % [7] human wealth
        'd'       % [8] dividend
        'r_w'     % [9] wealth return (detrended)
        'r'       % [10] equity return (detrended)
        'r_h'     % [11] human wealth return (detrend)
        'r_f'     % [12] risk-free return (non-detrended)
        'q'       % [13] marginal q
        'upsilon' % [14] time preference shock
        'nu'      % [15] investment shock
        'z'       % [16] technology shock
        'E_c'     % [17] expected c
        'E_r_w'   % [18] expected r_w
        'E_r'     % [19] expected r
        'd_c'     % [20] lag c
        'd_i'     % [21] lag i
        'd_w'     % [22] lag w
        'd_r_f'   % [23] lag r_f
        'sigma2'  % [24] volatility
        };

%%
% This section defines names for the structural shocks in the model. In
% the absence of measurment errors, the number of shocks cannot be less
% than the number of observable variables.

mshock = {
          'eps_u'  % [1] time preference
          'eps_n'  % [2] investment
          'eps_z'  % [3] technology
          'eps_s'  % [4] volatility
          };

%%
% This section defines names for the forecast errors in the model. These
% errors correspond to the expectation variables in the model.

mfore = {
         'c'    % [1] expected consumption
         'r_w'  % [2] expected aggregate wealth return
         'r'    % [3] expected equity return
         };

%%
% This section defines names for the observable variables used to estimate
% the model.

dvar = {
        'EXR'  % [1] log equity market excess return
        'RFR'  % [2] log real risk-free return
        'CGR'  % [3] real per capita consumption log growth rate
        'IGR'  % [4] real per capita investment log growth rate
        'WGR'  % [5] real per capita compensation log growth rate
        };

%% END