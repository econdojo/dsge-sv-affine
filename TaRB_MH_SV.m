function TaRB_MH_SV(M, burn, varargin)

% TaRB_MH_SV.m
%
% Implement TaRB-MH MCMC algorithm; see Chib & Ramamurthy (2010), 'TaRB
% MCMC methods with application to DSGE models', Journal of Econometrics
%
% Input
%
% M         sample size, including burn-ins (15000)
% burn      number of burn-ins  (5000)
% varargin  (optional) string/value pair
%           'pdof'   - proposal t degrees of freedom (15); pdof > 2
%           'prob'   - blocking probability (0.85)
%           'rec'    - number of iterations before new tailoring (1)
%           'nodisp' - suppress display if nonzero (1)
%
% Output
%
% none
%
% Written by Dave Rapach & Fei Tan, Saint Louis University
%
% Updated: 11-Jul-2018

clc
addpath('User', 'Toolbox', 'MacMex', 'WinMex');

%%
% This section initializes optional inputs and checks the validity of the
% inputs.

pdof   = 15;
prob   = 0.85;
rec    = 1;
nodisp = 1;

narg = length(varargin);

for k = 1:2:narg

    switch varargin{k}

        case 'pdof', pdof = varargin{ k+1 };
        case 'prob', prob = varargin{ k+1 };
        case 'rec', rec = varargin{ k+1 };
        case 'nodisp', nodisp = varargin{ k+1 };

        otherwise
            error('Unrecognized optional argument')

    end

end

if M-burn <= 500

    error('Number of draws after burn-in <= 500')

elseif pdof <= 2

    error('Student-t degrees of freedom <= 2')

elseif prob < 0 || prob > 1

    error('Blocking probability < 0 or > 1')

elseif rec < 0 || rec > 5

    error('Number of iterations before new tailoring < 0 or > 5')

end

%%
% This section loads parameters/variables/data and initializes entities
% for the MCMC algorithm.

[P, V] = ParVar;  % load parameters/variables

Y = importdata('User/Data.txt');  % import data

chain_theta = zeros(M, P.npara);      % Markov chain - parameters
chain_s2    = zeros(M, length(Y)+1);  % Markov chain - SV

invH = 1e-2*eye(P.npara);  % initial (-)inverse Hessian - proposal density

chain_theta( 1 , : ) = P.para( : , 1 )';  % start chains at initial values

s2 = ones(1, length(Y)+1)*P.para( P.sigma , 1 )^2/100;
s2 = PostSig(s2, P, V, Y, pdof);

chain_s2( 1 , : ) = s2;

pk_last = -PostKer1(P.para( : , 1 )', 1:P.npara, s2, P, V, Y);

rej1 = 0;  % number of overall rejections
rej2 = 0;  % number of recycle rejections
nb1  = 0;  % number of overall blocks
nb2  = 0;  % number of recycle blocks
next = 2;  % next round of randomization

%%
% This section implements the MCMC algorithm to simulate posterior draws.

tic;
progressbar('TaRB-MH MCMC in Progress')

for iter = 2:M

    fprintf('Iteration = %d\n', iter);
    fprintf('\n');

    % We compute randomized blocks for the parameters and determine
    % the number of sweeps for recycling the randomized blocks.

    if iter == next

        order = randperm(P.npara);  % random permutation

        k1 = 1;  % beginning index

        for k2 = 2:P.npara+1

            if rand <= prob && k2 <= P.npara

                continue

            end

            k1 = cat(2, k1, k2);  % next block

        end

        nb = length(k1) - 1;  % number of blocks w/in sweep

        theta_rec = cell(nb, 1);  % recycle theta
        Hess_rec  = cell(nb, 1);  % recycle (-)inverse Hessian

        n    = ceil(2*rec*rand+1e-3);  % recycle for next n sweeps
        next = iter + n;

    else

        nb2 = nb2 + nb;

        fprintf('Recycle\n');

    end

    % [1] The first step samples randomized blocks of parameters.
    
    for j = 1:nb

        block = order(k1( j ):k1( j+1 )-1);

        % When we don't recycle, we update by tailoring the proposal
        % density for each new randomized block using csminwel.m in
        % conjuction with PostKer1.m.

        if iter == next-n

            [~, theta_rec{ j }, ~, Hess_rec{ j }, ~, ~, rc] = ...
                csminwel('PostKer1', chain_theta( iter-1 , block ), ...
                invH( block , block ), [], 1e-5, 100, nodisp, block, ...
                s2, P, V, Y);

            if rc == 7

                Hess_rec{ j } = invH( block , block );

            end

            fprintf('Block %d csminwel return code = %d\n', j, rc);

        end
        
        % We use a Metropolis-Hastings step to draw the block of
        % parameters. We draw from the proposal density using mvt_rnd.m
        % and compute the move probability using MoveProb.m.

        theta = mvt_rnd_mex(theta_rec{ j }, Hess_rec{ j }, pdof, 1);

        [~, ~, alpha] = MoveProb('PostKer1', ...
            chain_theta( iter-1 , block ), theta, theta_rec{ j }, ...
            Hess_rec{ j }, pdof, pk_last, [], block, s2, P, V, Y);

        if rand > alpha

            chain_theta( iter , block ) = chain_theta( iter-1 , block );

            rej1 = rej1 + 1;

            if iter ~= next-n

                rej2 = rej2 + 1;

            end

        else

            chain_theta( iter , block ) = theta;

        end

        P.para( block , 1 ) = chain_theta( iter , block )';

    end

    % [2] The second step samples the stochastic variances.

    s2 = PostSig(s2, P, V, Y, pdof);

    chain_s2( iter , : ) = s2;

    pk_last = -PostKer1(P.para( : , 1 )', 1:P.npara, s2, P, V, Y);

    % We update and print sweep information.

    nb1 = nb1 + nb;

    progressbar(iter/M)

    fprintf('\n');
    disp('Parameter vector');
    disp(chain_theta( iter , : ));

    fprintf('Maximum variance = %f\n', max(chain_s2( iter , : )));
    fprintf('\n');
    fprintf('Recursive mean for psi (incl burn-ins) = %f\n', ...
        mean(chain_theta( 1:iter , 1 )));
    fprintf('Recursive sd for psi (incl burn-ins)   = %f\n', ...
        std(chain_theta( 1:iter , 1)));
    fprintf('\n');

    if iter >= burn+1

        fprintf('Recusrive mean for psi (post burn-ins) = %f\n', ...
            mean(chain_theta( burn+1:iter , 1 )));
        fprintf('Recursive sd for psi (post burn-ins)   = %f\n', ...
            std(chain_theta( burn+1:iter , 1)));
        fprintf('\n');

    end

end

time = datestr(toc/(24*60*60), 'DD:HH:MM:SS');

fprintf('Elapsed time = %s [dd:hh:mm:ss]\n', time);
fprintf('\n');

%%
% This section computes posterior means and 90% credible intervals for
% the parameters (after discarding the burn-in draws). The results are
% saved in the matrix stat. DiagPlot.m generates plots for convergence
% diagnostics.

chain_theta = chain_theta( burn+1:end , : );

chain_s2 = chain_s2( burn+1:end , : );  % need for TaRB_MH_SV.mat

M = M - burn;

stat = zeros(P.npara, 4);

stat( : , 1 )   = mean(chain_theta)';
stat( : , 2:3 ) = ProbBand(chain_theta, 0.9, 1);

disp('Para             Mean          90% Interval          Ineff')

for k = 1:P.npara

    stat( k , 4 ) = IneffFactor(chain_theta( : , k ), 200);

    fprintf('%s          %.3f         [%.3f, %.3f]        %.1f\n', ...
        P.name{ k }, stat( k , 1 ), stat( k , 2 ), stat( k , 3 ), ...
        stat( k , 4 ));
end

fprintf('\n');
fprintf('Overall rejection rate  =  %.1f %%\n', rej1/nb1*100);
fprintf('Recycle rejection rate  =  %.1f %%\n', rej2/nb2*100);
fprintf('Average inefficiency factor  =  %.1f\n', ...
    sum(stat( : , 4 ))/P.npara);
fprintf('Average number of blocks  =  %.1f\n', nb1/(M+burn-1));

save('User/TaRB_MH_SV.mat', 'chain_theta', 'stat', 'chain_s2');

DiagPlot(P, chain_theta)

rmpath('User', 'Toolbox', 'MacMex', 'WinMex');

%% END