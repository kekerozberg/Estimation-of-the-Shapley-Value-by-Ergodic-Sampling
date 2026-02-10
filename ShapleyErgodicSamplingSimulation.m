%% -----------------------------------------------------------------------
% ShapleyErgodicSamplingSimulation m-file simulates different estimates of
% the Shapley value of given games by ergodic and stratified sampling
% method
%
% For the result see Illes and Kerenyi [2019] Table 3
% ------------------------------------------------------------------------


%% -----------------------------------------------------------------------
% Programmed and implemented: Ferenc Illes
%                             Assistant Professor
%                             ferenc.illes@uni-corvinus.hu
%
%                             Peter Kerenyi
%                             Assistant Professor
%                             peter.kerenyi@uni-corvinus.hu
%
%                             Institute of Finance
%                             Corvinus Business School
%                             Corvinus University of Budapest
% ------------------------------------------------------------------------
% Paper:
% 1) Illes, F and Kerenyi, P (2019): Estimation of the Shapley Value by
% Ergodic Sampling
% ------------------------------------------------------------------------


%% Clear
% clear;

%% Parameters
% ------------------------------------------------------------------------

% % size of benchmark sample
m = 100000;
% 
% % number of sample for GreedyK2 optimalization algorithm
m1 = 50;
% 
% % random seed
% 
seed = 774;
% 
% 
% % number of players
n = 100;
% 
% % id of player
i = 1;
% 
% 
% % simulation number
simulationNumber = 4;
% 
% 
% %game = 'nonsymmetricvotinggame';       % Game 1
game = 'symmetricvotinggame';          % Game 2
% %game = 'shoesgame';                    % Game 3
% %game = 'airportgame';                  % Game 4
% game = 'minimumspanningtreegame';      % Game 5
% %game = 'bankruptcygame';               % Game 6
% %game = 'liabilitygame';                % Game 7
% %game = 'pairgame';                     % Game 8

rng(seed);

%% Simulation
% -------------------------------------------------------------------------

% size of sample
m2 = floor((m - (1/6) * m1 * n^2) ./ 2);

EstShapleyValues = nan(simulationNumber,1);

stratifiedErgodicEstShapleyValue = nan(simulationNumber,1);
stratifiedErgodicTotalComputationTime = nan(simulationNumber,1);
stratifiedErgodicElapsedTimeDetailedGreedy = nan(simulationNumber,1);
stratifiedErgodicElapsedTimeDetailedOptimalStratification = nan(simulationNumber,1);
stratifiedErgodicElapsedTimeDetailedStratifiedErgodicSampling = nan(simulationNumber,1);
elapsedTimeStratifiedErgodicSampling = nan(simulationNumber,1);


ergodicEstShapleyValue = nan(simulationNumber,1);
conditionalErgodicEstShapleyValue = nan(simulationNumber,1);
correlationGreedy = nan(simulationNumber,1);
ergodicTotalComputationTime = nan(simulationNumber,1);
ergodicElapsedTimeDetailedGreedy = nan(simulationNumber,1);
ergodicElapsedTimeDetailedErgodicSampling = nan(simulationNumber,1);
elapsedTimeErgodicSampling = nan(simulationNumber,1);

stratifiedEstShapleyValue = nan(simulationNumber,1);
conditionalStratifiedEstShapleyValue = nan(simulationNumber,1);
stratifiedTotalComputationTime = nan(simulationNumber,1);
stratifiedElapsedTimeDetailedOptimalStratification = nan(simulationNumber,1);
stratifiedElapsedTimeDetailedStratifiedSampling = nan(simulationNumber,1);
elapsedTimeStratifiedSampling = nan(simulationNumber,1);

randomEstShapleyValue = nan(simulationNumber,1);
conditionalRandomEstShapleyValue = nan(simulationNumber,1);
randomTotalComputationTime = nan(simulationNumber,1); 
randomElapsedTimeDetailedRandomSampling = nan(simulationNumber,1);
elapsedTimeRandomSampling = nan(simulationNumber,1);



parfor j=1:simulationNumber

    s = RandStream('CombRecursive','Seed',seed + j);
    RandStream.setGlobalStream(s);

    
    mainTimer = tic;
    [stratifiedErgodicEstShapleyValue(j,1), stratifiedErgodicTotalComputationTime(j,1), stratifiedErgodicElapsedTimeDetailedGreedy(j,1), stratifiedErgodicElapsedTimeDetailedOptimalStratification(j,1), stratifiedErgodicElapsedTimeDetailedStratifiedErgodicSampling(j,1)] = StratifiedErgodicShapley(m, m1, n, i, game);
    elapsedTimeStratifiedErgodicSampling(j,1) = toc(mainTimer);

    mainTimer = tic;
    [ergodicEstShapleyValue(j,1), conditionalErgodicEstShapleyValue(j,1), correlationGreedy(j,1), ergodicTotalComputationTime(j,1), ergodicElapsedTimeDetailedGreedy(j,1), ergodicElapsedTimeDetailedErgodicSampling(j,1)] = ConditionalErgodicShapley(2 .* m2, m1, n, i, game);
    elapsedTimeErgodicSampling(j,1) = toc(mainTimer);
    
    mainTimer = tic;
    [stratifiedEstShapleyValue(j,1), conditionalStratifiedEstShapleyValue(j,1), stratifiedTotalComputationTime(j,1), stratifiedElapsedTimeDetailedOptimalStratification(j,1), stratifiedElapsedTimeDetailedStratifiedSampling(j,1)] = StratifiedShapley(m, n, i, game);
    elapsedTimeStratifiedSampling(j,1) = toc(mainTimer);

    mainTimer = tic;
    [randomEstShapleyValue(j,1), conditionalRandomEstShapleyValue(j,1), randomTotalComputationTime(j,1), randomElapsedTimeDetailedRandomSampling(j,1)] = RandomShapley(m, n, i, game);
    elapsedTimeRandomSampling(j,1) = toc(mainTimer);

end

workspaceName = strcat("workspaces\_",num2str(simulationNumber),"_",game,"_",num2str(i),"_m_",num2str(m),"_m1_",num2str(m1),"___",replace(string(datetime),":","_"));
save(workspaceName);