% Author:    Abir Das and Anirban Chakraborty
% Date:      2014/11/25 15:23
% Revision:  1.0
% Copyright: Abir Das, Anirban Chakraborty and Amit K. RoyChowdhury, 2014

clear;
close all;
clc;

%% Parameter initialization and addition of paths
addpath('Pairwise_Similarity');
addpath('Results');
addpath('Results/GeneralizedNCR');

% Dataset name. This will be used to load pairwise wise similarity scores
% and to store intermediate variable and results. For example the pairwise
% similarity score will read from a file named -
% ['Pairwise_sim_Less_People_' dataset_name '.mat']; The optimization
% variables (intermediate) will be stored as ['OptVar_' dataset_name
% '.mat'] etc.
dataset_name = 'RAiD_6a'; % Use 'RAiD_6a' for generating the data to plot Fig 6(a), whereas use 'RAiD_6b' for generating the data to plot Fig 6(b)
% Number of trials
testset = 1:3;
% Camera numbers for this dataset
cameras = [1 2 3 4];
% Number of persons present in each camera. For ease of coding, in each
% test (consisting of random test people IDs), the last 4 persons are
% ommited from the last two cameras (in the example below)
numPersons = [20, 20, 12, 12]; % Use [20, 20, 16, 16] for generating the data to plot Fig 6(a), whereas use [20, 20, 12, 12] for generating the data to plot Fig 6(b)
% kk is the different 'k' values to try
kk = 0:0.05:0.4;

%% Load pairwise similarity scores and generate optimization cost function and constraints

% Load the pairwise similarity matrices. If the number of cameras be C,
% then this mat file will contain a cell array (named 'pairwise_sim') of
% size CxC, where diagonal as well as the lower diagonal elements of
% 'pairwise_sim' will be empty matrices. Other elements of 'pairwise_sim'
% contains the pairwise similarity scores between the persons in the
% corresponding cameras. As a concrete example, if the number of cameras be
% 4 (e.g. cam 1, cam 2, cam 3 and cam 4) and the number of persons be 20,
% 20, 15 and 15 respectively in these cameras with the number of tests or
% trials being 5, then 'pairwise_sim{1,2}' contains the 
% pairwise similarity scores between 20 persons between camera 1 and camera
% 2 for all the 5 tests. It is stored in a matrix of size 20x20x5, where
% the third dimension corresponds to the number of tests. Each 20x20 matrix
% is the pairwise similarity scores between the 20 persons for different
% testsets. Similarly 'pairwise_sim{1,3}' and 'pairwise_sim{2,3}' contains
% 20x15x5 matrices each containing pairwise similarity scores between
% camera 1-3 and camera 2-3.
pairwise_sim_filename = ['Pairwise_sim_Less_People_' dataset_name '.mat'];
load(fullfile('Pairwise_Similarity',pairwise_sim_filename))

% Get the different camera pairs
CPairs = combnk(cameras,2);
numCPairs = size(combnk(cameras,2), 1);
% Number of cameras
numCameras = length(cameras);

vectorSize_perpair = zeros(numCPairs, 1);
for i = 1:numCPairs
    n1 = CPairs(i,1);
    n2 = CPairs(i,2);
    vectorSize_perpair(i) = numPersons(n1)*numPersons(n2);
end

%% For different k values
for kCount = 1:length(kk)
    
    % Compute the cost function
    f = cell(length(testset),1);
    for iTSCount=testset
        for i = 1:numCPairs
            n1 = CPairs(i,1);
            n2 = CPairs(i,2);

            cc =  pairwise_sim{n1,n2}(1:numPersons(n1),1:numPersons(n2),iTSCount);
            f{iTSCount} = [f{iTSCount}; cc(:)];
        end;
        f{iTSCount} = kk(kCount).*ones(length(f{iTSCount}),1) - f{iTSCount};
    end
    disp('Cost vector done.');

    % Build the association/equality constraints
    A_assoc = [];
    vectorSize = size(f{1},1);
    for c = 1:numCPairs
        n1 = CPairs(c,1);
        n2 = CPairs(c,2); 
        % Compute the upper half of the matrix as given in equation (7) of the
        % supplementary. Each 'i' fills a row.
        for i = 1:numPersons(n1)
            v = zeros(1, vectorSize);
            for j = 1:numPersons(n2)
                v(1, sum(vectorSize_perpair(1:c-1))+(j-1)*numPersons(n1)+i) = 1;
            end
            A_assoc = [A_assoc; v];
        end
        % Compute the lower half of the matrix as given in equation (7) of the
        % supplementary. Each 'j' fills a row.
        for j = 1:numPersons(n2)
            v = zeros(1, vectorSize);
            for i = 1:numPersons(n1)
                v(1, sum(vectorSize_perpair(1:c-1))+(j-1)*numPersons(n1)+i) = 1;
            end
            A_assoc = [A_assoc; v];
        end
    end
    b_assoc = ones(size(A_assoc,1),1);
    disp('Association constraint matrix done.');

    % Get all possible triplets that can be formed out of the camera numbers.
    % Remember that this is not all possible triplets that can be formed by the
    % labels (x in paper), this is just all possible triplets that can be
    % formed by the camera numbers
    triplets = [];
    vectorSize_pertriplet = [];
    for r = 1:size(CPairs,1)
        cp = CPairs(r,1);
        cq = CPairs(r,2);
        others = setdiff(cameras,[cp,cq]);
        all_perms_tri = combntns(others, 1);
        for cr = all_perms_tri'
            triplets = [triplets; [cp,cr,cq]];
            vectorSize_pertriplet = [vectorSize_pertriplet; ...
                [numPersons(cp)*numPersons(cq)*numPersons(cr)]];
        end;
    end;

    % Fill up the loop/inequality constraints
    A_loop = zeros(sum(vectorSize_pertriplet), vectorSize);
    b_loop = zeros(size(A_loop,1), 1);
    for t = 1:size(triplets,1)
        % Camera p
        cp = triplets(t,1);
        % Camera r
        cr = triplets(t,2);
        % Camera q
        cq = triplets(t,3);
        % In 'CPairs' find the row number which contains the pair p-r
        cpr = find((CPairs(:,1)==cp & CPairs(:,2)==cr) | (CPairs(:,1)==cr & CPairs(:,2)==cp));
        % In 'CPairs' find the row number which contains the pair r-q
        crq = find((CPairs(:,1)==cr & CPairs(:,2)==cq) | (CPairs(:,1)==cq & CPairs(:,2)==cr));
        % In 'CPairs' find the row number which contains the pair p-q
        cpq = find((CPairs(:,1)==cp & CPairs(:,2)==cq) | (CPairs(:,1)==cq & CPairs(:,2)==cp));
        for i = 1:numPersons(cp)
            for k = 1:numPersons(cq)
                for j = 1:numPersons(cr)
                    v = zeros(1, vectorSize);
                    if cp < cr
                        v(1, sum(vectorSize_perpair(1:cpr-1))+(j-1)*numPersons(cp)+i) = 1;
                    else
                        v(1, sum(vectorSize_perpair(1:cpr-1))+(i-1)*numPersons(cr)+j) = 1;
                    end;
                    if cr < cq
                        v(1, sum(vectorSize_perpair(1:crq-1))+(k-1)*numPersons(cr)+j) = 1;
                    else
                        v(1, sum(vectorSize_perpair(1:crq-1))+(j-1)*numPersons(cq)+k) = 1;
                    end;
                    if cp < cq
                        v(1, sum(vectorSize_perpair(1:cpq-1))+(k-1)*numPersons(cp)+i) = -1;
                    else
                        v(1, sum(vectorSize_perpair(1:cpq-1))+(i-1)*numPersons(cq)+k) = -1;
                    end;
                    A_loop(sum(vectorSize_pertriplet(1:t-1)) + (i-1)*numPersons(cq)*numPersons(cr) + (k-1)*numPersons(cr) + j, :) = v;
                    b_loop(sum(vectorSize_pertriplet(1:t-1)) + (i-1)*numPersons(cq)*numPersons(cr) + (k-1)*numPersons(cr) + j, :) = 1;
                end;
            end;
        end;
    end;
    disp('Loop constraint matrix done.');

    % Now all the constraints are inequality, so club them
    A = [A_assoc; A_loop];
    b = [b_assoc; b_loop];

%% Run integer program to solve for the labels

    ReshapingIndices = [0; cumsum(vectorSize_perpair)];
    percentageAccuracy = zeros(size(CPairs,1),length(testset));
    fprintf('=====================\n');
    for iTSCount=testset
        t = tic;
        % Call cplex ip solver
        x = cplexbilp(f{iTSCount},A,b,[],[]);
        fprintf('\nkk = %0.2f Testset %d Optimization done in %.2f(s)\n', kk(kCount), iTSCount, toc(t));
        
        for i=1:size(CPairs,1)
            c1 = CPairs(i,1);
            c2 = CPairs(i,2);
            Allx{c1,c2}(:,:,iTSCount) = reshape(x(ReshapingIndices(i)+1:ReshapingIndices(i+1)),...
                numPersons(c1),numPersons(c2));
            % Calculate accuracy. Note that the assumption told in line
            % 26-28 (for numPersons) is employed here to avoid coding
            % complexity.
            % Calculate the number of correct matches
            correctMatch{c1,c2} = sum(diag(Allx{c1,c2}(:,:,iTSCount)));
            % Calculate the number of correct 'no matches'
            CorrectNoMatch{c1,c2} = 0;
            if size(Allx{c1,c2}(:,:,iTSCount),1) == size(Allx{c1,c2}(:,:,iTSCount),2)
                %
            elseif size(Allx{c1,c2}(:,:,iTSCount),1) > size(Allx{c1,c2}(:,:,iTSCount),2)
                RowsCorrToExtraPeople = Allx{c1,c2}(numPersons(c2)+1:numPersons(c1),:,iTSCount);
                CorrectNoMatch{c1,c2} = (numPersons(c1)-numPersons(c2))-length(find(RowsCorrToExtraPeople));
            else
                disp('Tell me');
            end
            
            correctResult = correctMatch{c1,c2} + CorrectNoMatch{c1,c2};
            percentageAccuracy(i,iTSCount) = correctResult/numPersons(c1)*100;
            fprintf('Camera %d-%d gave true positive %d ; true negative %d. Accuracy = %2.2f%%\n',...
                    c1,c2,correctMatch{c1,c2},CorrectNoMatch{c1,c2},percentageAccuracy(i,iTSCount))
        end
    end
    fprintf('\n');
    
    % Average Accuracy
    AvgPercentageAccuracy = mean(percentageAccuracy,2);
    for i=1:size(CPairs,1)
        c1 = CPairs(i,1);
        c2 = CPairs(i,2);
        fprintf('Camera %d-%d gave average accuracy = %2.2f%%\n',...
                    c1,c2,AvgPercentageAccuracy(i,1))
    end

    fprintf('\n\n');
    
end