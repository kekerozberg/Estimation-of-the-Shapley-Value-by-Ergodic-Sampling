function marginalContribution = getMarginalContribution(permutation, i, game)


switch lower(game)
    case 'nonsymmetricvotinggame'
        
        weightPlayers = [45, 41, 27, 26, 26, 25, 21, 17, 17, 14, 13, 13, ...
           12, 12, 12, 11, 10, 10, 10, 10, 9, 9, 9, 9, 8, 8, 7, 7, 7, 7, ...
           6, 6, 6, 6, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3,3];
        
        if size(permutation,2) == size(weightPlayers,2)
            marginalContribution = getNonSymmetricVotingGame(permutation, weightPlayers, i);
        end
        
        
    case 'symmetricvotinggame'
        
        marginalContribution = getSymmetricVotingGame(permutation, i);
        
        
    case 'shoesgame'
        
        marginalContribution = getShoesGame(permutation, i);
        
        
    case 'airportgame'
        
        weightPlayers = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, ...
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, ...
            4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, ...
            6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, ...
            8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, ...
            10, 10, 10, 10, 10];
     
        if size(permutation,2) == size(weightPlayers,2)
            marginalContribution = getAirportGame(permutation, weightPlayers, i);
        end
        
        
    case 'minimumspanningtreegame'

        marginalContribution = getMinimumSpanningTreeGame(permutation, i);

      
    case 'bankruptcygame'
        
        %% 6a
        weightPlayers = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, ...
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, ...
            4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, ...
            6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, ...
            8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, ...
            10, 10, 10, 10, 10];

        A = 200;
        
        %% 6b
%         weightPlayers = [9*ones(1,50), -8*ones(1,50)];
%         
%         A = 50;
        
     
        if size(permutation,2) == size(weightPlayers,2)
            marginalContribution = getBankruptcyGame(permutation, weightPlayers, A, i);
        end
        
        
    case 'liabilitygame'
        
        weightPlayers = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, ...
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, ...
            4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, ...
            6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, ...
            8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, ...
            10, 10, 10, 10, 10, 0];
        
        A = 200;
     
        if size(permutation,2) == size(weightPlayers,2)
            marginalContribution = getLiabilityGame(permutation, weightPlayers, A, i);
        end
        
    case 'pairgame'
        
        marginalContribution = getPairGame(permutation, i);
    
        
    otherwise
        marginalContribution = rand;
end
        
 

%% -----------------------------------------------------------------------
% getNonSymmetricVotingGame function computes the marginal contribution of 
% a given player for a given order in the Non symmetric voting game
%
%
% Syntax: marginalContribution = getNonSymmetricVotingGame(permutation, ...
%                                weightPlayers, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         weightPlayers - Weights of players (vector, number)
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
%
% ------------------------------------------------------------------------
function marginalContribution = getNonSymmetricVotingGame(permutation, weightPlayers, i)

sumAllPlayers = sum(weightPlayers);

orderPlayer = find(permutation == i);

sumPrecedingPlayers = sum(weightPlayers(permutation(1:(orderPlayer-1))));
sumPlayer = sum(weightPlayers(permutation(1:orderPlayer)));

if sumPlayer > (sumAllPlayers/2) && sumPrecedingPlayers <= (sumAllPlayers/2)
    marginalContribution = 1;
else
    marginalContribution = 0;
end



%% -----------------------------------------------------------------------
% getSymmetricVotingGame function computes the marginal contribution of a 
% given player for a given order in the Symmetric voting game
%
%
% Syntax: marginalContribution = getSymmetricVotingGame(permutation, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getSymmetricVotingGame(permutation, i)

numberPlayers = size(permutation,2);
orderPlayer = find(permutation == i);

if orderPlayer == floor((numberPlayers/2) + 1)
    marginalContribution = 1;
else
    marginalContribution = 0;
end



%% -----------------------------------------------------------------------
% getShoesGame function computes the marginal contribution of a given 
% player for a given order in the Shoes game
%
%
% Syntax: marginalContribution = getShoesGame(permutation, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         weightPlayers - Weights of players
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getShoesGame(permutation, i)

numberPlayers = size(permutation,2);
orderPlayer = find(permutation == i);
numberLeft = size(find(permutation(1:(orderPlayer-1)) <= floor(numberPlayers/2)), 2);
numberRight = size(find(permutation(1:(orderPlayer-1)) >= floor(numberPlayers/2+1)), 2);

if i <= (numberPlayers/2)
    marginalContribution = min(numberLeft+1, numberRight) - min(numberLeft, numberRight);
else
    marginalContribution = min(numberLeft, numberRight+1) - min(numberLeft, numberRight);
end



%% -----------------------------------------------------------------------
% getAirportGame function computes the marginal contribution of 
% a given player for a given order in the Airport game
%
%
% Syntax: marginalContribution = getAirportGame(permutation, ...
%                                weightPlayers, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         weightPlayers - Weights of players
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getAirportGame(permutation, weightPlayers, i)

orderPlayer = find(permutation == i);

if orderPlayer == 1
    marginalContribution = weightPlayers(permutation(1));
else
    marginalContribution = max(weightPlayers(permutation(1:orderPlayer))) - max(weightPlayers(permutation(1:(orderPlayer-1))));
end



%% -----------------------------------------------------------------------
% getMinimumSpanningTreeGame function computes the marginal contribution of 
% a given player for a given order in the Minimum spanning tree game
%
%
% Syntax: marginalContribution = getMinimumSpanningTreeGame(permutation, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getMinimumSpanningTreeGame(permutation, i)

numberPlayers = size(permutation,2);
orderPlayer = find(permutation == i);

if orderPlayer == numberPlayers
    marginalContribution = 1;
else
    leftNeighbor = i - 1;
    if leftNeighbor == 0
        leftNeighbor = numberPlayers;
    end

    rightNeighbor = i + 1;
    if rightNeighbor == (numberPlayers + 1)
        rightNeighbor = 1;
    end

    isLeft = any(permutation(1:1:(orderPlayer-1)) == leftNeighbor);
    isRight = any(permutation(1:1:(orderPlayer-1)) == rightNeighbor);

    if isLeft && isRight
        marginalContribution = -99;
    elseif isLeft || isRight
        marginalContribution = 1;
    else
        marginalContribution = 101;
    end
end



%% -----------------------------------------------------------------------
% getBankruptcyGame function computes the marginal contribution of 
% a given player for a given order in the Bankruptcy game
%
%
% Syntax: marginalContribution = getBankruptcyGame(permutation, ...
%                                weightPlayers, A, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         weightPlayers - Weights of players
%         A             - 
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getBankruptcyGame(permutation, weightPlayers, A, i)

orderPlayer = find(permutation == i);

if orderPlayer == size(permutation, 2)
    marginalContribution = max(0, A) - max(0, A - sum(weightPlayers(permutation(orderPlayer:end))));
else
    marginalContribution = max(0, A - sum(weightPlayers(permutation((orderPlayer+1):end)))) - max(0, A - sum(weightPlayers(permutation(orderPlayer:end)))); 
end



%% -----------------------------------------------------------------------
% getSquareGame function computes the marginal contribution of 
% a given player for a given order in the liability game
%
%
% Syntax: marginalContribution = getSquareGame(permutation, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         weightPlayers - Weights of players
%         A             - 
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getLiabilityGame(permutation, weightPlayers, A, i)

firmId = size(weightPlayers,2);
orderPlayer = find(permutation == i);
orderFirm = find(permutation == firmId);

if orderPlayer-1 < orderFirm
    v0 = max(0,A-sum(weightPlayers(permutation(orderPlayer:end))));
else
    v0 = min(A, sum(weightPlayers(permutation(1:(orderPlayer-1)))));
end

if orderPlayer < orderFirm
    v1 = max(0,A-sum(weightPlayers(permutation((orderPlayer+1):end))));
else
    v1 = min(A, sum(weightPlayers(permutation(1:orderPlayer))));
end

marginalContribution = v1-v0;



%% -----------------------------------------------------------------------
% getPairGame function computes the marginal contribution of a given player
% for a given order in the Pair game
%
%
% Syntax: marginalContribution = getPairGame(permutation, i)
%
% Inputs: permutation   - Order of players (vector, integer)
%         i             - ID of player (integer)
%
% Outputs: marginalContribution     - Marginal contribution of player
%                                     idPlayer in the order permutation
% ------------------------------------------------------------------------
function marginalContribution = getPairGame(permutation, i)

orderPlayer = find(permutation == i);
if mod(i,2) == 0
    pairPlayer = i-1;
else
    pairPlayer = i+1;
end
orderPairPlayer = find(permutation == pairPlayer);

if orderPlayer > orderPairPlayer
    marginalContribution =  1;
else
    marginalContribution =  0;
end