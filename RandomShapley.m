function [randomEstShapleyValue, conditionalRandomEstShapleyValue, totalComputationTime, elapsedTimeDetailedRandomSampling] = RandomShapley(m, n, i, game)

%% Variable declarations

Y = nan(m,1);
countOfOrder = zeros(1, n);
marginalContributionByOrder = zeros(1, n);
totalComputationTime = 0;


innerTimer = tic;
for j = 1:1:m
    randomPermutations = randperm(n);
    order = find(randomPermutations == i);
    Y(j,1) = getMarginalContribution(randomPermutations, i, game);
    totalComputationTime = totalComputationTime + 1;
    countOfOrder(order) = countOfOrder(order) + 1;
    marginalContributionByOrder(order) = marginalContributionByOrder(order) + Y(j, 1);
end

randomEstShapleyValue = mean(Y);
conditionalRandomEstShapleyValue = sum((marginalContributionByOrder ./ countOfOrder) .* (1/n));
elapsedTimeDetailedRandomSampling = toc(innerTimer);