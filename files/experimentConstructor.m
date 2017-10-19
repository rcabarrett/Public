%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 4a8b9c1r1v1R1V
%   - Output data structure decimation needed
%   - merge with presentedStims.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newstim, newcategory, sceneColours, newcolours, numFeatureLocations, numFeatures,featureFieldSize, catNum,featureAntiCorrelations,colourField,categoryFieldSize, catIndex, stimPos_x, stimPos_y, transfer, expName,LocationRelevance,LocationFeature] = experimentConstructor(structure, fieldModel, fitStruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimentConstructor %%
%
% Author: Jordan
% Date: Mar, 2013

%Brief description:
%Hard coded experiment parameters. ET3, SSHRC and FiveFour.

%Review: Huan, May 13...
%Verified: None.

% Inputs:
%
%       structure: integer. Value ranging from 1 to 3. Each value corresponds to a experiment.
%       fieldmodel: It has specific featureSpacing, categorySpacing and
%       halfSpace values.
%       fitStruct: specifies if fitting is in place.
%             fitStruct.fitting == 0 % Static: locations = funcRel
%             fitStruct.fitting == 1 % Randomized funcRels
%             fitStruct.fitting == 2 % Provided funcRels
%
%Outputs:
%       newstim: puts defined stimulus in random order and create 45 copies of them (21 for structure 3)
%       newcategory: puts defined category in random order and create 45 copies of them (21 for structure 3)
%       newcolour: puts defined colours in random order and create 45 copies of them (21 for structure 3)
%       sceneColours: extracts all 3 columns of colourField according to the value of the halfSpace and puts them into 3 separate rows.
%       numFeatures: total number of features with 2 features for each location
%       numFeatureLocations: number of locations for different structures. 3 for structure 1&2; 4 for structure 3.
%       featureFieldSize: The size of the feature field, calculated by the
%       product of the number of features and the feature spacing.
%       catNum: number of categories. 4 for structure 1,2; 2 for structure 3.
%       colourIndex: it is a number of features by 4 matrix.
%       colourField: 3 columns and featureFieldSize number of rows, each of which represents the dimension in the RGB colormap
%       categoryFieldSize: size of the category field and is calculated by
%       the number of categories multiplied by the category spacing
%       catIndex:
%       stimPos_x: x coordinate of the 3 locations with a radius of 17
%       stimPos_y: y coordinate of the 3 locations with a radius of 17
%       transfer: it puts pre-defined variables into new variables.
%       stimRadius: How far away each feature is from the origin.
%       expName: name of the experiment
%       LocationRelevance
%       LocationFeature
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if structure==1
    
    %% Constants
    expName = 'ET3';
    numberOfFitTrials = 360;
    numFeatureLocations = 3;
    numFeatures = numFeatureLocations * 2;
    catNum = 4;
    numStims = 8;
    blockNum = 360/numStims;
    stimRadius = 17;
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, 3, 5;
        1, 3, 6;
        1, 4, 5;
        1, 4, 6;
        2, 3, 5;
        2, 3, 6;
        2, 4, 5;
        2, 4, 6];
    
    category = [1;1;2;2;3;4;3;4];
    
    
    % The couples created by stimPos_x and stimPos_Y determine the
    % locations of the 3 features.
    
    stimPos_x = round(stimRadius*[ 0, cos(pi/6), -cos(pi/6) ] );
    stimPos_y = round(stimRadius*[ 1, -sin(pi/6), -sin(pi/6) ]);
    
    if fitStruct.fitting == 0
        newstim = [];
        newstim(:,fitStruct.funcRels(1)) = fitStruct.stimuli(1:numberOfFitTrials,2)+1;
        newstim(:,fitStruct.funcRels(2)) = fitStruct.stimuli(1:numberOfFitTrials,3)+3;
        newstim(:,fitStruct.funcRels(3)) = fitStruct.stimuli(1:numberOfFitTrials,4)+5;
        newcategory =  fitStruct.stimuli(1:numberOfFitTrials,1);
        newcolours = newstim;
        LocationRelevance.one=fitStruct.funcRels(1);LocationRelevance.two=fitStruct.funcRels(2);LocationRelevance.three=fitStruct.funcRels(3);
    else
        newcategory = [];
        newstim = [];
        newcolours = [];
        
        % This next loop takes the 8 stimuli that have been defined and created
        % 45 copies of them, randomly permuted such that you get an experiment
        % 360 trials in length, with stimuli appearing approximately 1 out of
        % every 8 trials but not in a particular set order. - HW
        
        if fitStruct.fitting == -1
            funcRels = [1 2 3];
            
        elseif fitStruct.fitting == 1
            funcRels = randperm(3);
        end
        
        stims = stims(:,funcRels);
        
        LocationRelevance.one=funcRels(1);
        LocationRelevance.two=funcRels(2);
        LocationRelevance.three=funcRels(3);
        
        for i = 1:blockNum
            randstim = randperm(numStims);
            newcategory = [newcategory;category(randstim',:)];
            newstim = [newstim;stims(randstim', :)];
            newcolours = [newcolours;stims(randstim', :)];
        end
        
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    
    featureAntiCorrelations =[];
    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1);
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        
        %These are what the columns of featureAntiCorrelations refer to:
        %[[feature value] [feature value location on the field] [anti-correlated feature value]
        %[anti-correlated feature value location on the field]
        %e.g.featureAntiCorrelations = [1 2 3 4 5 6;5 10 15 20 25 30;2 1 4 3 6 5; 10 5 20 15 30 25]';
        
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace];
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    transfer.switch=0;
    
end

%FiveFour Structure

if structure==3
    expName = '54';
    numFeatureLocations = 4;
    numFeatures = numFeatureLocations * 2;
    catNum = 2;
    numStims = 9;
    blockNum = 189/numStims;
    stimRadius = 17;
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    categoryFieldSize = catNum * categorySpacing;
    
    % In the next two matrices, the "0" is coded as larger numbers(2,4,6,8)
    % and the "1" is coded as smaller numbers(1,3,5,7) in feature values.
    % It is the opposite in struture 1.
    
    
    stims =[2, 4, 6, 7;
        2, 3, 6, 7;
        2, 3, 6, 8;
        2, 4, 5, 8;
        1, 4, 6, 8;
        2, 4, 5, 7;
        1, 4, 6, 7;
        1, 3, 5, 8;
        1, 3, 5, 7;];
    
    
    stimsTransfer =[2, 4, 6, 7;
        2, 3, 6, 7;
        2, 3, 6, 8;
        2, 4, 5, 8;
        1, 4, 6, 8;
        2, 4, 5, 7;
        1, 4, 6, 7;
        1, 3, 5, 8;
        1, 3, 5, 7;
        2, 3, 5, 8;
        2, 3, 5, 7;
        2, 4, 6, 8;
        1, 3, 6, 7;
        1, 4, 5, 8;
        1, 3, 6, 8;
        1, 4, 5, 7;];
    
    category = [1;1;1;1;1;2;2;2;2]; %F1 = Loc1, F2 = Loc2, F3 = Loc3
    categoryTransfer = [1;1;1;1;1;2;2;2;2;7;7;7;7;7;7;7];
    
    
    
    % The couples created by stimPos_x and stimPos_Y determine the
    % locations of the 4 features.
    
    stimPos_x = round(stimRadius*[ -cos(pi/4), cos(pi/4), cos(pi/4), -cos(pi/4) ] );
    stimPos_y = round(stimRadius*[ sin(pi/4) , sin(pi/4), -sin(pi/4), -sin(pi/4) ]);
    
    if fitStruct.fitting == 0
        newstim = [];
        newstim(:,fitStruct.funcRels(1)) = fitStruct.stimuli(1:numberOfFitTrials,2)+1;
        newstim(:,fitStruct.funcRels(2)) = fitStruct.stimuli(1:numberOfFitTrials,3)+3;
        newstim(:,fitStruct.funcRels(3)) = fitStruct.stimuli(1:numberOfFitTrials,4)+5;
        newcategory =  fitStruct.stimuli(1:numberOfFitTrials,1);
        newcolours = newstim;
        LocationRelevance.one=fitStruct.funcRels(1);
        LocationRelevance.two=fitStruct.funcRels(2);
        LocationRelevance.three=fitStruct.funcRels(3);
        LocationRelevance.four=fitStruct.funcRels(4);
    else
        newcategory = [];
        newstim = [];
        newcolours = [];
        
        % This next loop takes the 8 stimuli that have been defined and created
        % 45 copies of them, randomly permuted such that you get an experiment
        % 360 trials in length, with stimuli appearing approximately 1 out of
        % every 8 trials but not in a particular set order. - HW
        
        
        if fitStruct.fitting == -1
            funcRels = [1 2 3 4];
        elseif fitStruct.fitting == 1
            funcRels = randperm(4);
        end
        
        stims = stims(:,funcRels);
        LocationRelevance.one=funcRels(1);
        LocationRelevance.two=funcRels(2);
        LocationRelevance.three=funcRels(3);
        LocationRelevance.four=funcRels(4);
        
        for i = 1:blockNum
            randstim = randperm(numStims);
            newcategory = [newcategory;category(randstim',:)];
            newstim = [newstim;stims(randstim', :)];
            newcolours = [newcolours;stims(randstim', :)];
        end
        
        
    end
    
    transferCategories = [];
    transferStims = [];
    transferColours = [];
    
    % this loop is similar to the previous one except for using the transfer
    % stimuli.
    
    
    
    for i = 1:2
        transferCategories = [transferCategories;categoryTransfer];
        transferStims = [transferStims;stimsTransfer];
        transferColours = [transferColours;stimsTransfer];
    end
    
    stimsTransfer = stimsTransfer(:,funcRels);
    transferColours = stimsTransfer
    
    %still need the transfer stims for fitting.
    
    
    featureAntiCorrelations =[];
    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1); % it selects the values of first row and 1-8 columns in the colourField matrix.
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace]; % this equation creats a row of numbers that for odd i, it returns (i i i+1 i+1); for even i, it returns (i i i-1 i-1). Then the pattern repeats for through 1 o8 times to creat the featureAntiCorrelations matrix.
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    transfer.categories = transferCategories;
    transfer.stims = transferStims;
    transfer.colours = transferColours;
    transfer.switch = 1;
    
    
end


if structure==2
    expName = 'SSHRC';
    numFeatureLocations = 3;
    numFeatures = numFeatureLocations * 2;
    catNum = 4;
    numStims = 8;
    blockNum = 360/numStims;
    stimRadius = 17;
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize)
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, 3, 5;
        1, 3, 6;
        1, 4, 5;
        1, 4, 6;
        2, 3, 5;
        2, 3, 6;
        2, 4, 5;
        2, 4, 6];
    category = [1;1;2;2;3;3;4;4];
    
    
    % The couples created by stimPos_x and stimPos_Y determine the
    % locations of the 3 features.
    
    stimPos_x = round(stimRadius*[ 0, cos(pi/6), -cos(pi/6) ] );
    stimPos_y = round(stimRadius*[ 1, -sin(pi/6), -sin(pi/6) ]);
    
    if fitStruct.fitting == 0
        newstim = [];
        newstim(:,fitStruct.funcRels(1)) = fitStruct.stimuli(1:numberOfFitTrials,2)+1;
        newstim(:,fitStruct.funcRels(2)) = fitStruct.stimuli(1:numberOfFitTrials,3)+3;
        newstim(:,fitStruct.funcRels(3)) = fitStruct.stimuli(1:numberOfFitTrials,4)+5;
        newcategory =  fitStruct.stimuli(1:numberOfFitTrials,1);
        newcolours = newstim;
        LocationRelevance.one=fitStruct.funcRels(1);LocationRelevance.two=fitStruct.funcRels(2);LocationRelevance.three=fitStruct.funcRels(3);
    else
        newcategory = [];
        newstim = [];
        newcolours = [];
        
        % This next loop takes the 8 stimuli that have been defined and created
        % 45 copies of them, randomly permuted such that you get an experiment
        % 360 trials in length, with stimuli appearing approximately 1 out of
        % every 8 trials but not in a particular set order. - HW
        
        
        if fitStruct.fitting == -1
            funcRels = [1 2 3];
            
        elseif fitStruct.fitting == 1
            funcRels = randperm(3);
        end
        
        stims = stims(:,funcRels);
        
        LocationRelevance.one=funcRels(1);
        LocationRelevance.two=funcRels(2);
        LocationRelevance.three=funcRels(3);
        
        
        for i = 1:blockNum
            randstim = randperm(numStims);
            newcategory = [newcategory;category(randstim',:)];
            newstim = [newstim;stims(randstim', :)];
            newcolours = [newcolours;stims(randstim', :)];
        end
        
        
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    
    featureAntiCorrelations =[];
    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1);
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        
        %These are what the columns of featureAntiCorrelations refer to:
        %[[feature value] [feature value location on the field] [anti-correlated feature value]
        %[anti-correlated feature value location on the field]
        %e.g.featureAntiCorrelations = [1 2 3 4 5 6;5 10 15 20 25 30;2 1 4 3 6 5; 10 5 20 15 30 25]';
        
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace];
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    transfer.switch=0;
end

if structure==4
    expName = 'SHJ4';
    numFeatureLocations = 3;
    numFeatures = numFeatureLocations * 2;
    catNum = 2;
    numStims = 8;
    blockNum = 360/numStims;
    stimRadius = 17;
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize)
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    
    stims =[1, 3, 5;
        1, 4, 6;
        1, 4, 5;
        2, 4, 5;
        2, 4, 6;
        2, 3, 6;
        2, 3, 5;
        1, 3, 6];
    category = [1;1;1;1;2;2;2;2];
    
    
    if fitStruct.fitting == 0
        funcRels = [1 2 3];
        
    elseif fitStruct.fitting == 1
        funcRels = randperm(3);
        
    end
    stims = stims(:,funcRels);
    
    LocationRelevance.one=funcRels(1);
    LocationRelevance.two=funcRels(2);
    LocationRelevance.three=funcRels(3);
    
    
    
    % The couples created by stimPos_x and stimPos_Y determine the
    % locations of the 3 features.
    
    stimPos_x = round(stimRadius*[ 0, cos(pi/6), -cos(pi/6) ] );
    stimPos_y = round(stimRadius*[ 1, -sin(pi/6), -sin(pi/6) ]);
    
    if fitStruct.fitting == 3
        newstim = [];
        newstim(:,fitStruct.funcRels(1)) = fitStruct.stimuli(1:numberOfFitTrials,2)+1;
        newstim(:,fitStruct.funcRels(2)) = fitStruct.stimuli(1:numberOfFitTrials,3)+3;
        newstim(:,fitStruct.funcRels(3)) = fitStruct.stimuli(1:numberOfFitTrials,4)+5;
        newcategory =  fitStruct.stimuli(1:numberOfFitTrials,1);
        newcolours = newstim;
        LocationRelevance.one=fitStruct.funcRels(1);LocationRelevance.two=fitStruct.funcRels(2);LocationRelevance.three=fitStruct.funcRels(3);
    else
        newcategory = [];
        newstim = [];
        newcolours = [];
        
        % This next loop takes the 8 stimuli that have been defined and created
        % 45 copies of them, randomly permuted such that you get an experiment
        % 360 trials in length, with stimuli appearing approximately 1 out of
        % every 8 trials but not in a particular set order. - HW
        
        for i = 1:blockNum
            randstim = randperm(numStims);
            newcategory = [newcategory;category(randstim',:)];
            newstim = [newstim;stims(randstim', :)];
            newcolours = [newcolours;stims(randstim', :)];
        end
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    
    featureAntiCorrelations =[];
    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1);
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        
        %These are what the columns of featureAntiCorrelations refer to:
        %[[feature value] [feature value location on the field] [anti-correlated feature value]
        %[anti-correlated feature value location on the field]
        %e.g.featureAntiCorrelations = [1 2 3 4 5 6;5 10 15 20 25 30;2 1 4 3 6 5; 10 5 20 15 30 25]';
        
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace];
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    transfer.switch=0;
end

if structure==5
    
    
    
    %% Constants
    expName = 'CL';
    conditionCode = 1;
    numberOfFitTrials = 280;
    numFeatureLocations = 3;
    numFeatures = numFeatureLocations * 2;
    catNum = 4;
    numStims = 8;
    blockNum = 280/numStims;
    stimRadius = 17;
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    
    % Adapting from the CL experiment
    
    Feature1 = [0 0 0 0 1 1 1 1]';
    Feature2 = [0 0 1 1 0 0 1 1]';
    Feature3 = [0 1 0 1 0 1 0 1]';
    
    if conditionCode == 1 || conditionCode == 3 %Easy
        
        category = [1 1 2 2 3 3 4 4]';
        %categoryAlignment = ['A' 'A' 'B' 'B' 'C' 'C' 'D' 'D'];
        
    elseif conditionCode == 2 || conditionCode == 4 %Hard
        
        category = [1 2 3 4 5 6 7 8]';
        %categoryAlignment = [{'A1'} {'A2'} {'B1'} {'B2'} {'C1'} {'C2'} {'D1'} {'D2'}];
        
    end
    
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, 3, 5;
        1, 3, 6;
        1, 4, 5;
        1, 4, 6;
        2, 3, 5;
        2, 3, 6;
        2, 4, 5;
        2, 4, 6];
    
    
    
    % The couples created by stimPos_x and stimPos_Y determine the
    % locations of the 3 features.
    
    
    
    
    
    % Feature spacing condition:
    if conditionCode == 3 || conditionCode == 4 %Close
        stimRadius = 17; %400 % Distance of features from fixation
        stimRadius1 = 17; %400 %this puts the F1 feature at 21cm from the other features, as it is in the other condition.
    else % Far
        stimRadius1 = 17; %400
        stimRadius = 25.5; %400
    end
    
    featureSideLength = 3;
    fixationXY = [51 51];
    referenceAngle = 0;  % This makes reference to location 1. Angle = 0 is horizontally left from fixation, as is 360.
    
    
    Location1Radian = (referenceAngle+30/360)*pi*2;
    Location1Centre = [round(cos(Location1Radian)*stimRadius1)+fixationXY(1) (round(sin(Location1Radian)*stimRadius))+fixationXY(2)]; % contains X & Y of centre
    Location1Coords = [Location1Centre-(featureSideLength/2) Location1Centre+(featureSideLength/2)]'; % Given the centre, restore the X/Y coordinates for actual screen coords
    
    Location2Radian = (mod(referenceAngle+180, 360)/360)*pi*2; % 120 because there are 3 features, each 120 deg (of a circle) apart; mod 360 because we wrap if we're over
    Location2Centre = [(round(cos(Location2Radian)*stimRadius1))+fixationXY(1) (round(sin(Location2Radian)*stimRadius1))+fixationXY(2)];
    Location2Coords = [Location2Centre-(featureSideLength/2) Location2Centre+(featureSideLength/2)]';
    
    Location3Radian = (mod(referenceAngle+330, 360)/360)*pi*2;
    Location3Centre = [(round(cos(Location3Radian)*stimRadius1))+fixationXY(1) (round(sin(Location3Radian)*stimRadius))+fixationXY(2)];
    Location3Coords = [Location3Centre-(featureSideLength/2) Location3Centre+(featureSideLength/2)]';
    
    %AllLocationCoords = round([Location1Coords Location2Coords Location3Coords]); % indexing into column c gets you feature c
    AllLocationCoordsLAG-1 = round([Location1Centre;Location2Centre;Location3Centre]); % indexing into column c gets you feature c
    
    
    %     stimPos_x = round(stimRadius*[ 0, cos(pi/6), -cos(pi/6) ] );
    %     stimPos_y = round(stimRadius*[ 1, -sin(pi/6), -sin(pi/6) ]);
    
    
    stimPos_x = AllLocationCoordsLAG-1(:,1)';
    stimPos_y = AllLocationCoordsLAG-1(:,2)';
    
    
    
    
    
    
    
    
    
    if fitStruct.fitting == 0
        newstim = [];
        newstim(:,fitStruct.funcRels(1)) = fitStruct.stimuli(1:numberOfFitTrials,2)+1;
        newstim(:,fitStruct.funcRels(2)) = fitStruct.stimuli(1:numberOfFitTrials,3)+3;
        newstim(:,fitStruct.funcRels(3)) = fitStruct.stimuli(1:numberOfFitTrials,4)+5;
        newcategory =  fitStruct.stimuli(1:numberOfFitTrials,1);
        newcolours = newstim;
        LocationRelevance.one=fitStruct.funcRels(1);LocationRelevance.two=fitStruct.funcRels(2);LocationRelevance.three=fitStruct.funcRels(3);
    else
        newcategory = [];
        newstim = [];
        newcolours = [];
        
        % This next loop takes the 8 stimuli that have been defined and created
        % 45 copies of them, randomly permuted such that you get an experiment
        % 360 trials in length, with stimuli appearing approximately 1 out of
        % every 8 trials but not in a particular set order. - HW
        
        if fitStruct.fitting == -1
            funcRels = [1 2 3];
            
        elseif fitStruct.fitting == 1
            funcRels = randperm(3);
        end
        
        stims = stims(:,funcRels);
        
        LocationRelevance.one=funcRels(1);
        LocationRelevance.two=funcRels(2);
        LocationRelevance.three=funcRels(3);
        
        for i = 1:blockNum
            randstim = randperm(numStims);
            newcategory = [newcategory;category(randstim',:)];
            newstim = [newstim;stims(randstim', :)];
            newcolours = [newcolours;stims(randstim', :)];
        end
        
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    
    featureAntiCorrelations =[];
    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1);
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        
        %These are what the columns of featureAntiCorrelations refer to:
        %[[feature value] [feature value location on the field] [anti-correlated feature value]
        %[anti-correlated feature value location on the field]
        %e.g.featureAntiCorrelations = [1 2 3 4 5 6;5 10 15 20 25 30;2 1 4 3 6 5; 10 5 20 15 30 25]';
        
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace];
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    transfer.switch=0;
    
end


%% Based on Kruschke 2005

if structure==6
    
    expName = 'BlockingKruschke2005';
    numFeaturesTotal = 10;
    
    
    %% Early: Until 15/16 correct in two consecutive blocks, 3 block minimum
    
    numberOfFitTrials = 160; %10 block cap? check Kruschke2005
    numFeatureLocations = 2;
    numFeatures = 10;
    catNum = 4;
    numStims = 10;
    blockNum = 10;
    stimRadius = 17;
    newcategory = [];
    newstim = [];
    newcolours = [];
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    

    
    Feature1 = [0 -1 1 -1 2 -1 3 -1]';
    Feature2 = [-1 0 -1 1 -1 2 -1 3]';
%     Feature3 = [-1 -1 -1 -1 -1 -1 -1 -1]';
%     Feature4 = [-1 -1 -1 -1 -1 -1 -1 -1]';
    
    category = [1 1 2 2 3 3 4 4]';
    %categoryAlignment = ['X1' 'X1' 'X2' 'X2' 'Y1' 'Y1' 'Y2' 'Y2'];
    
    
    featureFieldSize = numFeaturesTotal * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, -1;
        -1, 1;
        2, -1;
        -1, 2;
        3, -1;
        -1, 3;
        4, -1;
        -1, 4];
    
    
    
    if numFeatureLocations ==4
    stimPos_x = round(stimRadius*[ -cos(pi/4), cos(pi/4), cos(pi/4), -cos(pi/4) ] );
    stimPos_y = round(stimRadius*[ sin(pi/4) , sin(pi/4), -sin(pi/4), -sin(pi/4) ]);
    
    else
    
        %2 feature location code. Commenting this out for now. July 17, 2017
    
        featureSideLength = 3;
        fixationXY = [0 0];

        Location1Centre = [-13 0]; % contains X & Y of centre
        Location1Coords = [Location1Centre-(featureSideLength/2) Location1Centre+(featureSideLength/2)]'; % Given the centre, restore the X/Y coordinates for actual screen coords

        Location2Centre = [13 0];
        Location2Coords = [Location2Centre-(featureSideLength/2) Location2Centre+(featureSideLength/2)]';

        AllLocationCoordsLAG-1 = round([Location1Centre;Location2Centre]); % indexing into column c gets you feature c


        stimPos_x = AllLocationCoordsLAG-1(:,1)';
        stimPos_y = AllLocationCoordsLAG-1(:,2)';

    end
    
    for i = 1:blockNum
        randstim = randperm(size(stims,1));
        newcategory = [newcategory;category(randstim',:)];
        newstim = [newstim;stims(randstim', :)];
        newcolours = [newcolours;stims(randstim', :)];
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    
    featureAntiCorrelations =[];
    for i = 1:numFeaturesTotal
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1);
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
        
        % There are no feature anti correlations for this structure because
        % all feature values will appear in all locations with any other
        % possible feature val...except itself...?
        featureAntiCorrelations = [featureAntiCorrelations; i (i*featureSpacing)-halfSpace i+(2*mod(i,2))-1 (featureSpacing*(i+(2*mod(i,2))-1))-halfSpace];
    end
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    
    %% Late: Until 15/16 correct in two consecutive blocks, 3 block minimum
    
    numberOfFitTrials = 160; %10 block cap? check Kruschke2005
    blockNum = 10;
    
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    Feature1 = [0 4 1 5 6 8 7 9]';
    Feature2 = [4 0 5 1 8 6 9 7]';
    
    category = [1 1 2 2 3 3 4 4]';
    %categoryAlignment = ['X1' 'X1' 'X2' 'X2' 'Y1' 'Y1' 'Y2' 'Y2'];
    
    
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, 5;
        5, 1;
        2, 6;
        6, 2;
        7, 9;
        9, 7;
        8, 10;
        10, 8];
    
    
    for i = 1:blockNum
        randstim = randperm(size(stims,1));
        newcategory = [newcategory;category(randstim',:)];
        newstim = [newstim;stims(randstim', :)];
        newcolours = [newcolours;stims(randstim', :)];
    end
    
    
    
    
    %% Test: 48 trials total
    
    
    numberOfFitTrials = 16+32; %2x8
    blockNum = 2;
    
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    
    %%
    
    Feature1 = [0 4 1 5 6 8 7 9]';
    Feature2 = [4 0 5 1 8 6 9 7]';
    
    category = [7 7 7 7 7 7 7 7]';
    

    
    stimsTransfer =[1, 5;
        5, 1;
        2, 6;
        6, 2;
        7, 9;
        9, 7;
        8, 10;
        10, 8];
    
    
    categoryTransfer = [7;7;7;7;7;7;7;7];
    
    transferCategories = [];
    transferStims = [];
    transferColours = [];
    
    
    for i = 1:2
        transferCategories = [transferCategories;categoryTransfer];
        transferStims = [transferStims;stimsTransfer];
        transferColours = [transferColours;stimsTransfer];
    end
    
    %A1 = 1, A2 = 2, F1 = 3, F2 = 4, B1 = 5, B2 = 6, C1 = 7, C2 = 8, D1 = 9, D2 = 10, 
    
    Feature1dba = [9 5 10 5 7 5 8 5]';
    Feature2dba = [5 9 5 10 5 7 5 8]';
    
    categorydba = [7 7 7 7 7 7 7 7]';
    
    Feature1dbb = [10 6 9 6 8 6 7 6]';
    Feature2dbb = [6 10 6 9 6 8 6 7]';
    
    categorydbb = [7 7 7 7 7 7 7 7]';
    
    Feature1aca = [1 7 1 8 1 9 1 10]';
    Feature2aca = [7 1 8 1 9 1 10 1]';
    
    categoryaca = [7 7 7 7 7 7 7 7]';
    
    Feature1acb = [2 8 2 7 2 10 2 9]';
    Feature2acb = [8 2 7 2 10 2 9 2]';
    
    categoryacb = [7 7 7 7 7 7 7 7]';
    
    
    Feature1 = [Feature1dba;Feature1dbb;Feature1aca;Feature1acb];
    Feature2 = [Feature2dba;Feature2dbb;Feature2aca;Feature2acb];
    
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    transferCategories = [transferCategories; categorydba;categorydbb;categoryaca;categoryacb];
    transferStims = [transferStims;[Feature1+1,Feature2+1]];
    transferColours = [transferColours;[Feature1+1,Feature2+1]];
    
    
    
    if fitStruct.fitting == 0
        funcRels = [1 2];
        
    elseif fitStruct.fitting == 1
        funcRels = randperm(2);
    else
        funcRels = [1 2];
        
    end
    
    newstim = newstim(:,funcRels);
    newcolours = newcolours(:,funcRels);
    transferStims = transferStims(:,funcRels);
    
    
    LocationRelevance.one=funcRels(1);
    LocationRelevance.two=funcRels(2);

    
    
    
    transfer.categories = transferCategories;
    transfer.stims = transferStims;
    transfer.colours = transferColours;
    transfer.switch = 1;
    
end





%% Simple Kamin blocking

if structure==7
    
    expName = 'SimpleBlocking';
    numFeaturesTotal = 4;
    
    
    %% Early: Until 15/16 correct in two consecutive blocks, 3 block minimum ? Kruschke2005 work
    
    numberOfFitTrials = 50; %10 block cap? check Kruschke2005
    numFeatureLocations = 4;
    numFeatures = 4;
    catNum = 2;
    numStims = 2;
    blockNum = 25;
    stimRadius = 17;
    newcategory = [];
    newstim = [];
    newcolours = [];
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    

    
    Feature1 = [0 -1]';
    Feature2 = [-1 -1]';
    Feature3 = [-1 2]';
    Feature4 = [-1 -1]';
    
    category = [1 2]';
    %categoryAlignment = ['A' 'B'];
    
    
    featureFieldSize = numFeaturesTotal * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, -1 ,-1, -1;
        -1, -1, 3, -1;];
    
    
    stimPos_x = round(stimRadius*[ -cos(pi/4), cos(pi/4), cos(pi/4), -cos(pi/4) ] );
    stimPos_y = round(stimRadius*[ sin(pi/4) , sin(pi/4), -sin(pi/4), -sin(pi/4) ]);
    

    
    for i = 1:blockNum
        randstim = randperm(size(stims,1));
        newcategory = [newcategory;category(randstim',:)];
        newstim = [newstim;stims(randstim', :)];
        newcolours = [newcolours;stims(randstim', :)];
    end
    
    
    
    %This loop gets the appropriate number of colours, equidistant in HSV
    %space.
    

    for i = 1:numFeatures
        sceneColours(1,i,1) = colourField((i*featureSpacing)-halfSpace,1); % it selects the values of first row and 1-8 columns in the colourField matrix.
        sceneColours(1,i,2) = colourField((i*featureSpacing)-halfSpace,2);
        sceneColours(1,i,3) = colourField((i*featureSpacing)-halfSpace,3);
    end
    featureAntiCorrelations = [1 1 3 3;2 2 4 4;3 3 1 1;4 4 2 2];
    
    
    catIndex = [];
    for i = 1:catNum
        catIndex = [catIndex (i*categorySpacing)-halfSpace];
    end
    
    
    
    %% Late
    
    numberOfFitTrials = 50; %10 block cap? check Kruschke2005
    blockNum = 25;
    
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    %%
    
    Feature1 = [0 -1]';
    Feature2 = [1 -1]';
    Feature3 = [-1 2]';
    Feature4 = [-1 3]';
    
    category = [1 2]';
    %categoryAlignment = ['A' 'B'];
    
    
    featureFieldSize = numFeaturesTotal * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    stims =[1, 2 ,-1, -1;
        -1, -1, 3, 4;];
    
    
    for i = 1:blockNum
        randstim = randperm(size(stims,1));
        newcategory = [newcategory;category(randstim',:)];
        newstim = [newstim;stims(randstim', :)];
        newcolours = [newcolours;stims(randstim', :)];
    end
    
    
    
    
    %% Test: 48 trials total
    
    
    numberOfFitTrials = 4; %2x8
    blockNum = 2;
    
    
    if fieldModel
        featureSpacing = 12;
        categorySpacing = 12;
        halfSpace = featureSpacing/2;
    else
        featureSpacing = 1;
        categorySpacing = 1;
        halfSpace = 0;
    end
    
    %%
    
    Feature1 = [-1 -1]';
    Feature2 = [1 -1]';
    Feature3 = [-1 -1]';
    Feature4 = [-1 3]';
    
    
    stimsTransfer =[1, 2 ,-1, -1;
        -1, -1, 3, 4];
    
    
    categoryTransfer = [7;7];
    
    transferCategories = [];
    transferStims = [];
    transferColours = [];
    
    
    for i = 1:2
        transferCategories = [transferCategories;categoryTransfer];
        transferStims = [transferStims;stimsTransfer];
        transferColours = [transferColours;stimsTransfer];
    end
    

    
    featureFieldSize = numFeatures * featureSpacing; %add 1 for padding
    colourField = hsv(featureFieldSize);
    
    categoryFieldSize = catNum * categorySpacing;
    
    
    
    if fitStruct.fitting == 0
        funcRels = [1 2 3 4];
        
    elseif fitStruct.fitting == 1
        funcRels = randperm(4);
    else
        funcRels = [1 2 3 4];
        
    end
    
    newstim = newstim(:,funcRels);
    newcolours = newcolours(:,funcRels);
    transferStims = transferStims(:,funcRels);
    
    
    LocationRelevance.one=funcRels(1);
    LocationRelevance.two=funcRels(2);
    LocationRelevance.three=funcRels(3);
    LocationRelevance.four=funcRels(4);


    
    transfer.categories = transferCategories;
    transfer.stims = transferStims;
    transfer.colours = transferColours;
    transfer.switch = 1;
    
end


featureOptions = {}; %Describes the feature colours used.
for i = 1:2:numFeatures
    [iscolor clrs] = fuzzycolor(hsv2rgb(squeeze(sceneColours(:,i,:))'));
    [iscolor1 clrs1] = fuzzycolor(hsv2rgb(squeeze(sceneColours(:,i+1,:))'));
    featureOptions{i} = strcat(clrs(find(iscolor==1,1)), '/', clrs1(find(iscolor1==1,1)));
    
end
%This is a SQL table thing. Would be nice to write the colours to SQL but
%currently this causes LAG-1 to break so just nominal numbers for now.


if numFeatureLocations ==3
    %Location1Feature=featureOptions{1};Location2Feature=featureOptions{3};Location3Feature=featureOptions{5};
    LocationFeature.one=1;LocationFeature.two=1;LocationFeature.three=1;
elseif numFeatureLocations ==4
    %Location1Feature=featureOptions{1};Location2Feature=featureOptions{3};Location3Feature=featureOptions{5};Location4Feature=featureOptions{7};
    LocationFeature.one=1;LocationFeature.two=1;LocationFeature.three=1;LocationFeature.four=1;
elseif numFeatureLocations ==2
    LocationFeature.one=1;LocationFeature.two=1;
    
end








end