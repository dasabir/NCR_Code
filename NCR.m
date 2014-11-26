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
addpath('Results/NCR');

% Dataset name. This will be used to load pairwise wise similarity scores
% and to store intermediate variable and results. For example the pairwise
% similarity score will read from a file named - ['Pairwise_sim_'
% dataset_name '.mat']; The optimization variables (intermediate) will be
% stored as ['OptVar_' dataset_name '.mat'] etc.
dataset_name = 'WARD';
% Number of trials
testset = 1:5;
% Camera numbers for this dataset
cameras = [1 2 3];
% Number of persons in each camera
numPersons = 35;

% Parameters specific to plots
labelFontSize = 14;
labelFontWeight = 'bold';
fullScreen = false;
%% Load pairwise similarity scores and generate optimization cost function and constraints

% Load the pairwise similarity matrices. If the number of cameras be C,
% then this mat file will contain a cell array (named 'pairwise_sim') of
% size CxC, where diagonal as well as the lower diagonal elements of
% 'pairwise_sim' will be empty matrices. Other elements of 'pairwise_sim'
% contains the pairwise similarity scores between the persons in the
% corresponding cameras. As a concrete example, if the number of cameras be
% 3 (e.g. cam 1, cam 2 and cam 3) and the number of persons be 35 and the
% number of tests or trials be 5, then 'pairwise_sim{1,2}' contains the
% pairwise similarity scores between 35 persons between camera 1 and camera
% 2 for all the 5 tests. It is stored in a matrix of size 35x35x5, where
% the third dimension corresponds to the number of tests. Each 35x35 matrix
% is the pairwise similarity scores between the 35 persons for different
% testsets. Similarly 'pairwise_sim{1,3}' and 'pairwise_sim{2,3}' contains
% 35x35x5 matrices each containing pairwise similarity scores between
% camera 1-3 and camera 2-3. All other elements of the cell array are empty
% matrices
pairwise_sim_filename = ['Pairwise_sim_' dataset_name '.mat'];
load(fullfile('Pairwise_Similarity',pairwise_sim_filename))


% Get the different camera pairs
CPairs = combnk(cameras,2);
numCPairs = size(combnk(cameras,2), 1);
% Number of cameras
numCameras = length(cameras);


% If saved matrices for this dataset is there load it, otherwise compute
% those
if exist(['Results\NCR\OptVar_' dataset_name '.mat'],'file')
    load(['Results\NCR\OptVar_' dataset_name '.mat'])
else
    % Build the association/equality constraints
    Aeq = [];
    for c = 1:numCPairs
        % Compute the upper half of the matrix as given in equation (7) of the
        % supplementary. Each 'i' fills a row.
        for i = 1:numPersons
            v = zeros(1, numCPairs*numPersons^2);
            for j = 1:numPersons
                v(1, (c-1)*numPersons^2+(j-1)*numPersons+i) = 1;
            end
            Aeq = [Aeq; v];
        end
        % Compute the lower half of the matrix as given in equation (7) of the
        % supplementary. Each 'j' fills a row.
        for j = 1:numPersons
            v = zeros(1, numCPairs*numPersons^2);
            for i = 1:numPersons
                v(1, (c-1)*numPersons^2+(j-1)*numPersons+i) = 1;
            end
            Aeq = [Aeq; v];
        end
    end
    beq = ones(size(Aeq,1),1);
    disp('Association constraint matrices done.');

    % Get all possible triplets that can be formed out of the camera numbers.
    % Remember that this is not all possible triplets that can be formed by the
    % labels (x in paper), this is just all possible triplets that can be
    % formed by the camera numbers
    triplets = [];
    for r = 1:size(CPairs,1)
        cp = CPairs(r,1);
        cq = CPairs(r,2);
        others = setdiff(cameras,[cp,cq]);
        all_perms_tri = combnk(others, 1);
        for cr = all_perms_tri'
            triplets = [triplets; [cp,cr,cq]];
        end;
    end;

    % Fill up the loop/inequality constraints
    B = zeros(size(triplets,1)*numPersons^3, numCPairs*numPersons^2);
    b = zeros(size(B,1), 1);
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
        % Fill up the values as in equation 12 of the supplementary
        for i = 1:numPersons
            for k = 1:numPersons
                for j = 1:numPersons
                    v = zeros(1, numCPairs*numPersons^2);
                    if cp < cr
                        v(1, (cpr-1)*numPersons^2+(j-1)*numPersons+i) = 1;
                    else
                        v(1, (cpr-1)*numPersons^2+(i-1)*numPersons+j) = 1;
                    end;
                    if cr < cq
                        v(1, (crq-1)*numPersons^2+(k-1)*numPersons+j) = 1;
                    else
                        v(1, (crq-1)*numPersons^2+(j-1)*numPersons+k) = 1;
                    end;
                    if cp < cq
                        v(1, (cpq-1)*numPersons^2+(k-1)*numPersons+i) = -1;
                    else
                        v(1, (cpq-1)*numPersons^2+(i-1)*numPersons+k) = -1;
                    end;
                    B((t-1)*numPersons^3 + (i-1)*numPersons^2 + (k-1)*numPersons + j, :) = v;
                    b((t-1)*numPersons^3 + (i-1)*numPersons^2 + (k-1)*numPersons + j, :) = 1;
                end;
            end;
        end;
    end;
    disp('Loop constraint matrices done.');

    % Compute the cost function
    f = cell(length(testset),1);
    for iTSCount=testset
        for i = 1:numCPairs
            n1 = CPairs(i,1);
            n2 = CPairs(i,2);
            cc =  pairwise_sim{n1,n2}(1:numPersons,1:numPersons,iTSCount);
            f{iTSCount} = [f{iTSCount}; cc(:)];
        end;
        f{iTSCount} = ones(length(f{iTSCount}),1) - f{iTSCount};
    end
    disp('Cost vectors done.');

    save(['Results\NCR\OptVar_' dataset_name '.mat'],'f','B','Aeq','b','beq','numPersons', 'numCameras', 'CPairs', 'testset', 'cameras', '-v7.3');
    %clear
end

    
%% Run integer program to solve for the labels
if exist(['Results\NCR\Res_' dataset_name '.mat'],'file')
    load(['Results\NCR\Res_' dataset_name '.mat']);
else
    for iTSCount=testset

        t = tic;
        % Call cplex ip solver
        x = cplexbilp(f{iTSCount},B,b,Aeq,beq);

        fprintf('Testset %d done in %.2f(s)\n', iTSCount, toc(t));

        for i=1:size(CPairs,1)
            c1 = CPairs(i,1);
            c2 = CPairs(i,2);
            Allx{c1,c2}(:,:,iTSCount) = reshape(x((i-1)*numPersons*numPersons+1:i*numPersons*numPersons),numPersons,numPersons);
            % Calculate accuracy
            correctMatch{c1,c2} = sum(diag(Allx{c1,c2}(:,:,iTSCount)));
            percentageAccuracy{c1,c2} = sum(diag(Allx{c1,c2}(:,:,iTSCount)))/numPersons*100;
            fprintf('Camera %d-%d gave %d out of %d correct matches. Accuracy = %2.2f%%\n',...
                c1,c2,correctMatch{c1,c2},numPersons,percentageAccuracy{c1,c2})
        end
    end

    % Save x's
    save(['Results\NCR\Res_' dataset_name '.mat'],'Allx','correctMatch','percentageAccuracy');
end

%% Create the new pairwise similarity score matrices and CMCs

if exist(['Results\NCR\ps_CMCs_new_' dataset_name '.mat'],'file')
 load(['Results\NCR\ps_CMCs_new_' dataset_name '.mat']);
else
    AllMeanCMCs_New = zeros(numPersons,size(CPairs,1));
    new_pairwise_sim = pairwise_sim;
    for i=1:size(CPairs,1)

        avgCMC = zeros(numPersons,1);
        c1 = CPairs(i,1);
        c2 = CPairs(i,2);
        for j=testset
            [row,col] = find(Allx{c1,c2}(:,:,j));
            for k=1:size(row,1)
                new_pairwise_sim{c1,c2}(row(k),col(k),j) = 1;
            end

            CMC = zeros(numPersons,1);
            avgMatchScoreOfThisTest = new_pairwise_sim{c1,c2}(:,:,j);
            for k = 1:numPersons
                sortedAvgMatchScore = sort(avgMatchScoreOfThisTest(k,:),'descend');
                true_order = find(sortedAvgMatchScore==avgMatchScoreOfThisTest(k,k));
                % The situation when 1 people in camera 1 has same score with
                % more than one people in the second camera. In that case take
                % the order which is least to break the tie. As the scores are
                % sorted, the orders will be sequential
                if length(true_order) > 1
                    true_order = true_order(1);
                end
                CMC(true_order) = CMC(true_order)+1;
            end
            CMC = CMC*100/length(CMC);
            CMC = cumsum(CMC);
            avgCMC = avgCMC + CMC;
        end
        AllMeanCMCs_New(:,i) = avgCMC/length(testset);
    end

    save(['Results\NCR\ps_CMCs_new_' dataset_name '.mat'],'new_pairwise_sim','AllMeanCMCs_New');
end


%% Generate the CMCs from the original pairwise similarity scores

if exist(['Results\NCR\CMCs_old_' dataset_name '.mat'],'file')
 load(['Results\NCR\CMCs_old_' dataset_name '.mat']);
else
    AllMeanCMCs_Old = zeros(numPersons,size(CPairs,1));
    for i=1:size(CPairs,1)

        avgCMC = zeros(numPersons,1);
        c1 = CPairs(i,1);
        c2 = CPairs(i,2);
        for j=testset
            CMC = zeros(numPersons,1);
            avgMatchScoreOfThisTest = pairwise_sim{c1,c2}(:,:,j);
            for k = 1:numPersons
                sortedAvgMatchScore = sort(avgMatchScoreOfThisTest(k,:),'descend');
                true_order = find(sortedAvgMatchScore==avgMatchScoreOfThisTest(k,k));
                % The situation when 1 people in camera 1 has same score with
                % more than one people in the second camera. In that case take
                % the order which is least to break the tie. As the scores are
                % sorted, the orders will be sequential
                if length(true_order) > 1
                    true_order = true_order(1);
                end
                CMC(true_order) = CMC(true_order)+1;
            end
            CMC = CMC*100/length(CMC);
            CMC = cumsum(CMC);
            avgCMC = avgCMC + CMC;
        end
        AllMeanCMCs_Old(:,i) = avgCMC/length(testset);
    end
    save(['Results\NCR\CMCs_old_' dataset_name '.mat'],'AllMeanCMCs_Old');
end

%% Plot the new and old CMCs

% Plot the CMCs for original pairwise similarity scores
for i=1:size(CPairs,1)    
    hCMC(i) = figure;
    hold on
    axis square
    grid on
    customLegend{1} = 'Original';
    
    plot(AllMeanCMCs_Old(:,i), 'Color', 'k', ...
                    'LineStyle', '-', ...
                    'LineWidth', 2, ...
                    'Marker', '>', ...
                    'MarkerSize', 4, ...
                    'MarkerEdgeColor', 'k');
end

% Plot CMCs for NCR given pairwise similarity scores
for i=1:size(CPairs,1)
    figure(hCMC(i));
    customLegend{2} = 'NCR';
    
    plot(AllMeanCMCs_New(:,i), 'Color', 'r', ...
                    'LineStyle', '-', ...
                    'LineWidth', 2, ...
                    'Marker', 'o', ...
                    'MarkerSize', 4, ...
                    'MarkerEdgeColor', 'r');
end

for i=1:size(CPairs,1)
    figure(hCMC(i));
    % Set plot title
    c1 = CPairs(i,1); c2 = CPairs(i,2);
    title({'Cumulative Matching Characteristic (CMC)'; ['Camera Pair ' ...
        num2str(c1) ' - ' num2str(c2)]});
    xlabel('Rank Score');
    ylabel('Recognition Percentage');
    xlim([1 numPersons]);
    ylim([1 100]);
    set(gca,'xtick',[0:5:numPersons])
    set(gca,'ytick',[0:10:100])
    %Add legend to plot
    legend(customLegend, 'Location', 'SouthEast');
    hold off;
end

% Set custom style and show
for i=1:size(CPairs,1)
    set_label_sytle(hCMC(i), labelFontSize, labelFontWeight, fullScreen);
end




