% *-----------------------------------*
% |   DSGE-SV-affine MATLAB Toolbox   |
% |                                   |
% |   (c) David E. Rapach & Fei Tan   |
% |   Last modified: 14-May-2019      |
% *-----------------------------------*
%
% I. GENERAL INFORMATION
%
% The MATLAB(R) program in this package implements an efficient Bayesian
% algorithm for estimating DSGE models with stochastic volatility. An
% affine solution procedure for the DSGE model allows for the estimation of
% risk premia. The Bayesian algorithm is described in our paper, "Estimating
% Risk Premia in DSGE Models with Stochastic Volatility." The toolbox is
% provided for researchers and comes with no performance guarantees.
% 
% II. CONTENTS
%
% The DSGE-SV-affine/ folder contains a set of .m programs:
%
%   See also
%   DiagPlot.m
%   MoveProb.m
%   ParVar.m
%   PostKer1.m
%   PostKer2.m
%   PostSig.m
%   Prior.m
%   TaRB_MH_SV.m
%   TestBench.m
%
% The DSGE-SV-affine folder also contains a set of subfolders:
%
%   User/     - data file and a set of user input .m programs
%   Toolbox/  - set of general-purpose .m programs
%   MacMex/   - C source MEX functions for Mac OS users
%   WinMex/   - C source MEX functions for Windows OS users
% 
% III. EXAMPLE
%
%   1. Set the MATLAB current working folder to DSGE-SV-affine/.
%
%   2. Modify Data.txt and the other files in the User/ subfolder as
%      needed.
%
%   3. At the MATLAB command line, call the TaRB_MH_SV.m function by
%      typing, for example,
%      >> TaRB_MH_SV(15000, 5000, 'pdof', 10, 'prob', 0.6, 'rec', 3)
%
%   4. The parameter and stochastic volatility draws (without the burn-
%      ins) are stored in the MATLAB data file, TaRB_DSGE_SV.mat, which
%      is saved to the User/ subfolder. All MEX functions needed can be
%      generated using the MATLAB Coder by calling, e.g., at the command
%      line:
%      >> TestBench
% 
% IV. CONTACT
%
% If you experience difficulty using this package or have comments for
% improvement, please feel free to contact us via e-mail:
%
% David E. Rapach            Fei Tan
% david.rapach@slu.edu       fei.tan@slu.edu
% 
% V. REVISION NOTES
%
% 1. 21-Jul-2018 update: original version
%
% 2. 14-May-2019 update: corrects some bugs found by Kadir Ozen
%
% END

clc
help ReadMe