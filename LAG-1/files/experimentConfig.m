% NOT CURRENTLY IN USE. READ NOTES IN MAIN FILES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 6w0r0v
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% experimentConfig.m %%
%
% Author: Jordan Barnes
% Date: 2016
%
% Description:
%   Configure LAG-1 simulation parameters
%
%  
%  This script requires eval or assign statements. Unused for now.
%    S = coder.load('experimentConfig.mat')
%    assigns= structfun(@(f) [f ' = S.' f '; '],S,'uniformoutput',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Toggle global noise:
noise = true;
%
%   Logging:
logging = 1;
%
%   Turn the GUI on: % deprecate because of line 56?
visualize = 0;
%
%	Use feature and category fields
fieldModel = 0;
%
%   Turns on GUI at a particular trial number. Use this ONLY IF visualize is set to 0
startVisualize = 1;
endVisualize = 2;


%% Additional experiment parameters

tau = 10; %Time scaling factor
deltaT = 2; %Characteristic time scale.
endFit = 0; %Fitting finished? Set to NaN if early model termination.
learning_rate = 0.000013% 0.000009;  1.3000e-05
softMax = 1; %Luce decision rule on or off.
antiHebbCorr = 1; %turns anti-Hebbian learning on for corrective.
trialNum = 1; % First Trial
gridRefinement = 0.25; %Factor allowing the model to scale spatially. i.e. you want things like stimuli and kernels to scale with field size.
spatialFieldSize = ceil(gridRefinement * 201); %Arbitrary initial scale.
visualFieldSize = 2*(spatialFieldSize-1)+1;
spatialHalfSize = (spatialFieldSize-1)/2; %allow you to call center of field
visualHalfSize = (visualFieldSize-1)/2;  % ditto
targetSize = 25*gridRefinement;  %stimulus width
badTimestepNum = 5000; % If it's this long model will abort.
badFixationNum = 30; % Max fixation count before fail.

%This is for testing purposes.
%rng(1);

%Debugging
profiling = 0; % For model performance optimization using the profiler.
maxFields = 1; % Display the maximum values for particular fields.

save experimentConfig.mat
