function [stratifiedErgodicEstShapleyValue, totalComputationTime, stratifiedErgodicElapsedTimeDetailedGreedy, stratifiedErgodicElapsedTimeDetailedOptimalStratification, stratifiedErgodicElapsedTimeDetailedStratifiedErgodicSampling] = StratifiedErgodicShapley(m, m1, n, i, game)

w = 0.5;


%% Variable declarations

elapsedTime = nan(3,1);
totalComputationTime = 0;
MC = zeros(n,2);
MC2 = zeros(n,2);

m_il_exp_const = floor((w * m)/(n));
marginalContributions = nan(n, m1 + m_il_exp_const);


%% Stage 1
% Variable declarations
innerTimer = tic;
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
    marginalContributions(positionsi(j), MC(positionsi(j),2)) = greedyX(j);
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

elapsedTime(1) = toc(innerTimer);

%% Define optimal stratification
% Generate the missing samples for std calculation
innerTimer = tic;

% Variable declarations
std_i = nan(n, 1);
usedSampleSizes = MC(:,2);
missingSampleSizes = max(0, m_il_exp_const - usedSampleSizes);
newSampleSizes = (w * m - totalComputationTime) .* (missingSampleSizes ./ sum(missingSampleSizes));

for l = 1:n   % position l
    for j=1:newSampleSizes(l)   % sample j
        randomPermutation = randperm(n);
        positioni = randomPermutation == i;
        playerOnOrderl = randomPermutation(l);
        randomPermutation(positioni) = playerOnOrderl;
        randomPermutation(l) = i;
        Y = getMarginalContribution(randomPermutation, i, game);
        MC(l,1) = MC(l,1) + Y;
        MC(l,2) = MC(l,2) + 1;
        MC2(l,1) = MC2(l,1) + Y;
        MC2(l,2) = MC2(l,2) + 1;
        marginalContributions(l,MC(l,2)) = Y;
        totalComputationTime = totalComputationTime + 1;
    end
    std_i(l) = std(marginalContributions(l,:),'omitmissing');
end
usedSampleSizes = MC(:,2);

% Estimate stds
m_i = m .* (std_i ./ sum(std_i));
missingSampleSizes = max(0, m_i - usedSampleSizes(l));
if sum(missingSampleSizes,'omitmissing') == 0
    missingSampleSizes(:) = m_il_exp_const;
end
newSampleSizes = (m - totalComputationTime) .* (missingSampleSizes ./ sum(missingSampleSizes));
elapsedTime(2) = toc(innerTimer);



%% Stage 2
%  Generate stratified and ergodic sample
innerTimer = tic;

mx = nan(n,3);
mx(:,1) = min(newSampleSizes, newSampleSizes(pairsPlayers));
mx(:,2) = newSampleSizes - mx(:,1);

for l=1:n
    p = pairsPlayers(l);
    for j=1:mx(l,1)
        randomPermutations1 = randperm(n);
        orderOfPlayer = randomPermutations1 == i;
        playerOnOrderl = randomPermutations1(l);
        randomPermutations1(orderOfPlayer) = playerOnOrderl;
        randomPermutations1(l) = i;
        tempMC = getMarginalContribution(randomPermutations1, i, game);
        totalComputationTime = totalComputationTime + 1;
        MC2(l,1) = MC2(l,1) + tempMC;
        MC2(l,2) = MC2(l,2) + 1;
        if l~=p
           randomPermutations2 = randomPermutations1(pairsPlayers);
           tempMC = getMarginalContribution(randomPermutations2, i, game);
           totalComputationTime = totalComputationTime + 1;
           MC2(p,1) = MC2(p,1) + tempMC;
           MC2(p,2) = MC2(p,2) + 1;
        end
    end
    mx(p,1) = mx(p,1) - mx(l,1);
    
    for j=1:mx(l,2)
        randomPermutations1 = randperm(n);
        orderOfPlayer = randomPermutations1 == i;
        playerOnOrderl = randomPermutations1(l);
        randomPermutations1(orderOfPlayer) = playerOnOrderl;
        randomPermutations1(l) = i;
        tempMC = getMarginalContribution(randomPermutations1, i, game);
        totalComputationTime = totalComputationTime + 1;
        MC2(l,1) = MC2(l,1) + tempMC;
        MC2(l,2) = MC2(l,2) + 1;
    end
end



Sh_i = MC2(:,1) ./ MC2(:,2);

stratifiedErgodicEstShapleyValue = (1/n) * sum(Sh_i,'omitmissing');
elapsedTime(3) = toc(innerTimer);

stratifiedErgodicElapsedTimeDetailedGreedy = elapsedTime(1);
stratifiedErgodicElapsedTimeDetailedOptimalStratification = elapsedTime(2);
stratifiedErgodicElapsedTimeDetailedStratifiedErgodicSampling = elapsedTime(3);
