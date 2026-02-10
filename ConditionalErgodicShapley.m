function [ergodicEstShapleyValue, conditionalErgodicEstShapleyValue, correlationGreedy, totalComputationTime, elapsedTimeDetailedGreedy, elapsedTimeDetailedErgodicSampling] = ConditionalErgodicShapley(m2, m1, n, i, game)

%% Variable declarations
% Order of f function
K = 2;

% ID of player
if ~exist('i','var')
    i = 1;
end

if ~exist('game','var')
    game = 'minimumspanningtreegame';
end


L = round(m2/K);

Y = nan(L,K);
countOfOrder = zeros(1, n);
marginalContributionByOrder = zeros(1, n);
totalComputationTime = 0;
MC = zeros(n,2);
MC2 = zeros(n,2);



% Stage 1 of ergodic method to find a suitable transformation
innerTimer = tic;

% Variable declarations
positionsi = nan(m1,1);
pairsPlayers = (1:1:n)';
randomPermutations = cell(m1,1);
tranPermutations = cell(m1,1);
greedyX = nan(m1,1);
greedyY = nan(m1,1);
BetaOrg = nan(n,n);

%% Greedy
% Generate independently and uniformly random orders of players
for j = 1:1:m1
    randomPermutations{j} = randperm(n);
    positionsi(j) = find(randomPermutations{j,1} == i);
    greedyX(j) = getMarginalContribution(randomPermutations{j,1}, i, game);
    totalComputationTime = totalComputationTime + 1;
    MC(positionsi(j),1) = MC(positionsi(j),1)+ greedyX(j);
    MC(positionsi(j),2) = MC(positionsi(j),2) + 1;
end


for r = 1:1:n       % position r
    for s = r:1:n       % position s
        for j = 1:1:m1      % random permuation j
            valuer = randomPermutations{j}(r);
            values = randomPermutations{j}(s);
            positioni = positionsi(j);
            if positioni <= s && positioni >= r
                tranPermutations{j} = randomPermutations{j};
                tranPermutations{j}(r) = values;
                tranPermutations{j}(s) = valuer;
                positioni = find(tranPermutations{j} == i);
                greedyY(j) = getMarginalContribution(tranPermutations{j}, i, game);
                totalComputationTime = totalComputationTime + 1;
            else
                greedyY(j) = greedyX(j);
            end
            MC(positioni,1) = MC(positioni,1) + greedyY(j);
            MC(positioni,2) = MC(positioni,2) + 1;
            covMatrix = cov(greedyX,greedyY);
            BetaOrg(r,s) = covMatrix(1,2); % Covariance of data sets
        end
    end
end

Beta = BetaOrg;

% Search edges with minimal weights
for j=1:1:(ceil(n/2))  
    [~, I] = min(Beta(:));
    [I_row, I_col] = ind2sub(size(Beta),I);
    pairsPlayers(I_row) = I_col;
    pairsPlayers(I_col) = I_row;

    Beta(I_row,:) = nan;
    Beta(I_col,:) = nan;
    Beta(:,I_row) = nan;
    Beta(:,I_col) = nan;
end

elapsedTimeDetailedGreedy = toc(innerTimer);


%% Stage 2
%  Generate ergodic sample

innerTimer = tic;
for j = 1:1:L
    randomPermutations1 = randperm(n);
    randomPermutations2 = randomPermutations1(pairsPlayers);
    order1 = find(randomPermutations1 == i);
    order2 = find(randomPermutations2 == i);
    
    countOfOrder(order1) = countOfOrder(order1) + 1;
    Y(j, 1) = getMarginalContribution(randomPermutations1, i, game);
    totalComputationTime = totalComputationTime + 1;
    MC2(order1,1) = MC2(order1,1) + Y(j, 1);
    MC2(order1,2) = MC2(order1,2) + 1;

    marginalContributionByOrder(order1) = marginalContributionByOrder(order1) + Y(j, 1);

    countOfOrder(order2) = countOfOrder(order2) + 1;
    Y(j, 2) = getMarginalContribution(randomPermutations2, i, game);
    totalComputationTime = totalComputationTime + 1;
    MC2(order2,1) = MC2(order2,1) + Y(j, 2);
    MC2(order2,2) = MC2(order2,2) + 1;
    marginalContributionByOrder(order2) = marginalContributionByOrder(order2) + Y(j, 2);    
end

Sh_i = MC2(:,1) ./ MC2(:,2);

ergodicEstShapleyValue = mean(mean(Y));
conditionalErgodicEstShapleyValue = (1/n) * sum(Sh_i,'omitmissing');
correlationGreedy = corr(Y(:,1), Y(:,2));
elapsedTimeDetailedErgodicSampling = toc(innerTimer);