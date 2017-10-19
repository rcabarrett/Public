%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 6b7c1r1v1R1V
%	- refactor using Cosivina
%	- line count decimation
%	- data structure parsimony
%	- equation parsing and labeling in latex
%   - version 2.0 fixes needed:
%       - fbButtonDetection = fbButtonDetection(:,:,7);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [subjectNumber, trialLogger, endFit] = TwoDSimulator(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAG-1
%
% Author: Jordan Barnes, Calen Walshe, Mark Blair, Paul Tupper
% Date: July 2012 - 2016
% Reviewed: October 2013. Mark Blair & Paul Tupper.
% Verified: Summer 2016. Mark Blar & Paul Tupper.
%
% 1. preconditions and postconditions.
% 2. model components
% 3. additional experiment parameters
% 4. Structure and analysis setup
% 5. Cog-Bio parameters
% 6. Field, neuron and weight initializations
% 7. feedbackType, category, parameter fit settings
% 8. run-time parameters
% 9. Visualization
% 10. Trial initialization
% 11. Trial start
% 12. Feedback
% 13. Hebbian learning
% 14. Saccade and rotation mechanics
% 15. Neural update calculations

%
%% 1. pre-conditions and postconditions
%
% preconditions
%
% 1. An experiment structure has been defined in experimentConstructor.m
% 2. The control switches have been set.
%
%
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

%
%   Category structure type
if ~isempty(varargin)
    if ischar(varargin{1})
        structure = str2num(varargin{1}); %Fitting model directly.
    else
        structure = varargin{1}; %Calling model from another script
    end
else
    structure =1; %1=5to1,3=5/4
end
%
%   Set the learning feedbackType:
%   feedbackType = 1 = corrective feedback
%   feedbackType = 2 = learning only on correct trials
%   feedbackType = 3 = anti-Hebbian learning on error trials
%   feedbackType = 4 = anti-Hebbian and Hebbian learning on error trials
if length(varargin) == 2
    feedbackType = varargin{2};
else
    feedbackType = 1;
end

%% 2. Model components
% field_v - retinal field, space times feature, high resolution
% field_a - spatial attention field
% field_s - saccade motor field
% neurons_f - feature neurons
% neurons_c - categorization neurons
% neuron_r - saccade reset neuron
% neuron_x - fixation neuron, excites foveal position in field_a
% neuron_g - gaze-change neuron, inhibits foveal position in field_a
% neuron_i - inhibition of return neuron, inhibits previous foveal position
% neuron_click - decision to click neuron
% neuron_imp - scalar exponential factor
%%


%% 3. Additional experiment parameters

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


%% 4. Structure and analysis setup

[s,id]=system('echo $PBS_JOBID');
[a b c d] = regexp(id,'\d')
a = strcat(cell2mat(d))
[time] = clock; %[year, month, day, hour, minute, second]
rand('twister',sum(1e5*clock))
stamp = [num2str(time(6)) num2str(time(5)) num2str(time(4)) num2str(time(3)) num2str(time(2)) num2str(time(1))];
stamp(find(stamp=='.'))='';
subjectNumber = num2str(stamp);
if length(id) > 1
    if length(varargin)>2 %if we've got an actual subject number we're fitting
        fitStruct = varargin{5};
        subjectNumber = str2num([num2str(fitStruct.subject) subjectNumber(1:4) a(end-3:end)])
    else
        subjectNumber = str2num([subjectNumber(1:4) a(6:end)])
    end
else
    subjectNumber = str2num(subjectNumber(1:10));
end

%Trial level logging
CorrectResponse = [];
resultsStimulus = [];
feedbackOff = 0; %Set to 1 for experiment transfer phases where there's no FB.
trainingCategories = []; %Stores the categories before being wiped by transfer categories.
recordSaccStart = []; %For quick analyses looking at fixation durations and saccade times.
recordSaccEnd = []; %For quick analyses looking at fixation durations and saccade times.
cpDetector = 0; % CP=criterion point. Flags once criterion point is met.
TETTime = 0; %total time on experiment (annoying eye tracking variable name "Tobii Eye Tracker")
TETTimeTimeSteps = 0;
accuracyLevels = {[0.1:1:25]'}; %labels the levels of gamma/temperature assessed if desired.


%% Get experiment structure


% A varargin of length 5 indicates that we're fitting LAG1. This is
% currently called by TwoDSimulator(structure,feedback, paramNames, paramVec,fitStruct)
if length(varargin) == 5
    fitStruct = varargin{5};
    optimizationModeFlag =1;
else
    fitStruct.fitting=-1; %Location == funcRels;
    optimizationModeFlag=0;
end

% Obtain the parameters of the experiment.
[stimulus, category, sceneColors, colours, locationNum, featureNum, featureFieldSize, catNum, featureAntiCorrelations, colourField, catFieldSize, catIndex, stimPos_x, stimPos_y,transfer,expName, LocationRelevance, LocationFeature] = experimentConstructor(structure,fieldModel,fitStruct);

maxTrials = length(stimulus);
visualLayers = featureNum + 1;
totalLocations = locationNum + 1;


%% Feedback button placement

%Visual field
stimPos_x = [stimPos_x -18]; %Appending the fb locations
stimPos_y = [stimPos_y 20];

%Attention field 

spatialFBposition = [-18 20];

fbButtonSize = 32*gridRefinement;


%% Formatting log display

trialLogger(1:21) = [{'Category neurons at decision time'},...
    {'Sigmoided Category neurons at decision time'}, {'Feature neurons at decision time'},...
    {'Sigmoided Feature neurons at decision time'}, {'Weights at decision time'},...
    {'Sigmoided weights at decision time'}, {'Phase 2 Reaction time'},...
    {'Phase 4 Reaction time'}, {'Raw weight changes'},...
    {'Phase 2 Fixations'}, {'Phase 2 Fixation durations'},...
    {'Phase 4 Fixations'}, {'Phase 4 Fixation durations'}, {'Correct Category'},...
    {'Response'}, {'Total Trial Time'}, {'Saccade Times'},{'Accuracy'}, {'Accuracy Levels'}, {'Prob Fix Irrel'}, {'Fixation Changes'}];

trialLogger(2,1:17) = {[]};
trialLogger(2,21) = {[]};
trialLogger(1,24) = {'Func_Rels'};
trialLogger(2,24) = {LocationRelevance};


%% Formatting for SQL import

if logging
    wtHist = [];
    gazeHeader = {'Subject', 'ID', 'TETTime', 'CursorX', 'CursorY', 'ValidityLeftEye', 'ValidityRightEye', 'TrialID', 'TrialPhase'};
    txt=sprintf('%s\t',gazeHeader{:});
    txt(end)='';
    dlmwrite(['./dft' expName 'GazeLvl-' num2str(subjectNumber) '-1.gazedata'],txt ,'');
    GazeLvl = [];
    if exist(['./dft' expName 'ExpLvl.txt'])==0
        
        if locationNum ==3
            expHeader = {'Subject', 'Cond', 'Location1Feature', 'Location2Feature', 'Location3Feature', 'Location1Relevance', 'Location2Relevance', 'Location3Relevance'};
        elseif locationNum ==4
            expHeader = {'Subject', 'Cond', 'Location1Feature', 'Location2Feature', 'Location3Feature', 'Location4Feature', 'Location1Relevance', 'Location2Relevance', 'Location3Relevance', 'Location4Relevance'};
        elseif locationNum==2
            expHeader = {'Subject', 'Cond', 'Location1Feature', 'Location2Feature', 'Location1Relevance', 'Location2Relevance'};
            
        end
        
        txt=sprintf('%s\t',expHeader{:});
        txt(end)='';
        dlmwrite(['./dft' expName 'ExpLvl.txt'],txt,'');
    end
    ExpLvl = [];
    if locationNum ==3
        trialHeader={'Subject', 'TrialID', 'Feature1Value', 'Feature2Value', 'Feature3Value', 'CorrectResponse', 'Response', 'TrialAccuracy',	'StimulusRT', 'FeedbackRT',	'StartTime', 'FixationOnset', 'StimulusOnset','ColourOnset', 'FixCrossLocationX', 'FixCrossLocationY', 'FeedbackOnset'};
    elseif  locationNum ==4
        trialHeader={'Subject', 'TrialID', 'Feature1Value', 'Feature2Value', 'Feature3Value', 'Feature4Value','CorrectResponse', 'Response', 'TrialAccuracy',	'StimulusRT', 'FeedbackRT',	'StartTime', 'FixationOnset', 'StimulusOnset','ColourOnset', 'FixCrossLocationX', 'FixCrossLocationY', 'FeedbackOnset'};
    elseif locationNum==2
        trialHeader={'Subject', 'TrialID', 'Feature1Value', 'Feature2Value', 'CorrectResponse', 'Response', 'TrialAccuracy',	'StimulusRT', 'FeedbackRT',	'StartTime', 'FixationOnset', 'StimulusOnset','ColourOnset', 'FixCrossLocationX', 'FixCrossLocationY', 'FeedbackOnset'};
    end
    txt=sprintf('%s\t',trialHeader{:});
    txt(end)='';
    dlmwrite(['./dft' expName 'TrialLvl-' num2str(subjectNumber) '.txt'],txt, '');
    
end

%%


%% 5. Cog-Bio parameters


%cog_bio_parameters; %could try this with some kind of struct explode?

% q scales all temporal noise levels
q = 0.1;

% resting levels "h_#", betas "b_#", and noise levels "q_#", offset of
% sigma "inflec_#" for all fields
% _v - retinal field, space times feature, high resolution
% _a - spatial attention field
% _s - saccade motor field
% _f - feature neurons
% _pref - feature expectation neurons
% _c - categorization neurons
% _r - saccade reset neuron
% _x - fixation neuron, excites foveal position in field_a
% _g - gaze-change neuron, inhibits foveal position in field_a
% _i - inhibition of return neuron, inhibits previous foveal position
% _wT - matrix of weights for learning connection between features and
% categories
% _click - decision to click neuron
% _imp - scalar exponential factor

h_v_1 = -4; beta_v_10 = 2; q_v_5 = 0.05 * q; inflec_v_11=2;
h_f_80 = -2; beta_f_85 = 1; q_f_95 = 0.2 * q; inflec_f_86 = 0; q_spatial_f_96 = 0.2;
h_pref_106 = -3; beta_pref_15 = 0.5; q_pref_111 = 1 * q; inflec_pref_16=0; q_spatial_pref_112 =1;
h_a_18 = -0.5;beta_a_32 = 1; q_a_39 = 1 * q; inflec_a_33 =0; q_spatial_a_40= 100;
h_s_47 = -2; beta_s_43 = 4; q_s_57 = 2 * q; inflec_s_44=0;q_spatial_s_58= 10;
h_r_75 = -5; beta_r_26 = 2; q_r_93 = 0.05 * q; inflec_r_27 = 0; q_spatial_r_94 = 1;
h_x_68 = -3.5; beta_x_50 = 0.5; q_x_91 = 0.1 * q; inflec_x_51 =0; q_spatial_x_92 = 1;
h_g_63 = -1; beta_g_30 = 3; q_g_89 = 0.01 * q; inflec_g_31=0; q_spatial_g_90 = 1;
h_c_97 = -0.5; beta_c_101 = 1.37; q_c_103 = 0.2 * q; inflec_c_102 = 0.15; q_spatial_c_104 = 0.2;
%h_i = -2; beta_i = 1; q_i = 1 * q; inflec_i=0;
h_wT = -0.86; beta_wT = 3; inflec_wT = 0.5;
h_imp_122 = 1; beta_imp_120 = .01; inflec_imp_121 = 300; q_imp_124 = 0.1 * q; q_spatial_imp_125 = 1;
h_click_107 = -5; beta_click_116 = 2.5; q_click_118 = 1 * q; inflec_click_117=1; q_spatial_click_119 = 0.1
h_f_80bButton = -5; beta_f_85bButton = 2.5; q_fbButton = q; inflec_f_86bButton=3;
h_f_80bButtonExp = 0.5; beta_f_85bButtonExp_23 = 1; q_fbButtonExp = q; inflec_f_86bButtonExp_24=0.75;


% c: strength of interaction kernel
% sigma: width of interaction kernel
% gi: strength of global inhibition
% exc = excitatory and inh = inhibitory

% NOTE: the field labels are reversed from what you might think
% e.g., c_fv = connection strength *FROM* visual field *TO* features.

sigma_q_spatial_17 = 20*gridRefinement; % For the noise kernel.

% Visual field parameters
c_vinp_1_2 = 4;
c_vr_1_3 = 3; %Gain on visual field suppression from saccade neuron.
c_vFBnovelty_1_4 = 1;
c_vv_exc_1_6 = 1;
sigma_vv_exc_spt_1_7 = 7*gridRefinement;
c_vv_inh_1_8 = 1;
sigma_vv_inh_s_47pt_1_9 = 15*gridRefinement;
c_va_gi_1_12 = 0.001;
c_vpref_1_13 = 3; %
c_vpref_gi_1_14 = 0.04;
sigma_vpref = 4*gridRefinement; %Only in field version
c_va_1_5 = 0; %Unused
sigma_va_1_4 = 5*gridRefinement; %Unused

% Attention Field parameters
c_apref_2_19 =0.2;
c_ax_2_20 = 5;
c_ai_2_21 = -0.3;
sigma_aIn_2_22=10*gridRefinement;
c_afbButton_2_25 = 5;
c_ag_2_28 = -5;
c_ar_2_29 = -10;
c_aa_exc_2_34 = 5;
sigma_aa_exc_2_35 = 20*gridRefinement;
c_aa_inh_2_36 = 8;%Increasing this allows you to see the convolution size.
sigma_aa_inh_2_37 = 30*gridRefinement;
c_aa_gi_2_38 = 0.2;
c_as_2_41 = 15;
sigma_as_2_42 = 10*gridRefinement;
c_av_2_45 = 70;
sigma_av_2_46 = 5*gridRefinement;
c_aFBnovelty = 1; %unused

% Saccade motor field parameters
c_sx_3_48 = -50;
sigma_ss_exc_3_49 = 20*gridRefinement;
c_sg_3_52 = 30;
c_sa_3_53 = 7;
sigma_sa_3_54 = 15*gridRefinement;
sigma_sfov_3_55 = 4; %This is actually used in calculating input to field_s
c_sfov_3_56 = 3;
c_ss_exc_3_59 = 5;
c_ss_inh_3_60 = 8;
sigma_ss_inh_3_61 = 30*gridRefinement;
c_ss_gi_3_3_62 = 5;
c_sa_gi_4 = 0.0; %unused
c_sr_gi_5 = 0; % unused

% Gaze change neuron parameters
c_gr_4_64 = -1;
c_gx_4_65 = 0.2;
c_gg_4_66 = 0.5;
c_gimp_4_67 = 0.0001;

% Fixation neuron parameters
c_xg_5_69 =-0.4;
c_xr_5_70 = -3;
c_xx_5_71 = 3;
c_xfeature_5_72 = 5;
beta_ft_x_73 = 1; % applied to foveal detection values
inflec_ft_x_74 = 1; % applied to foveal detection values

% Saccade Initiation/Reset Neuron parameters
c_rx_6_76 = -10;
c_rg_6_77 = 4;
c_rs_6_78 = 8;
c_rr_6_79 = 1.4;

% Feature neuron parameters
c_ff_inh_7_81 = 4;
c_fv_7_82 = 500;
c_ff_7_83 = 4.6;
c_ff_gi_7_84 = 0.1;
ft_det_fbButton_gi = c_ff_gi_7_84;
beta_ft_input_87 = 1; % feature detection sums get sigmoided.
inflec_ft_input_88 = 0; % feature detection sums get sigmoided.

% Category neuron parameters
c_cf_8_98 = 4.6;
c_cc_exc_8_99 =1.2; % Very unorthodox usage. The category is excited only by the other categories.
c_cc_inh_8_100 = 3.5;
c_cgain_8_105 = 1;

% Feature expectation neuron parameters
c_prefgain_9_105 = 0.2;
c_preff_9_108 = 0.01; % This can act as a kind of IOR.
c_prefc_9_109 = 2;
c_prefpref_exc_9_110 = 1.5;

% Click decision neuron parameters
c_clickc_10_113 = 4;
c_clickimp_10_114 = 1.6;
c_clickx_10_115 = -1.2;

% Impatience parameters
c_impTimer_11_123 = 0.001;
fixation_impatience_exponent = 1.7% 1.6;
trial_impatience_exponent = 1.75% 1.65;

% Feedback dynamics
c_xfbButton_changes = 1;


% Miscellaneous model parameters
temperature=0.07; %otherwise known as gamma. This sets how deterministic the model's responses should be.
decay_rate = 5; %higher decay constant here actually means lower decay (denominator)
antiHebbRate = 3;
CategoryChoiceThreshold = 0.8;
advanceTrialThreshold = 0.98;
c_xcategory_changes = 50;
c_xfeature_changes = 150;
c_inp_1 = 0.45;
fovealSize = ceil(20*gridRefinement);
fovealHalfSize = ceil(fovealSize/2);
sigma_fov_mask = 3;
feedback_c_strength = 0.3;
impatienceChangeStrength = -100;
iorTime = round(600*tau/deltaT); %timestep units. deltaT and tau are flipped because iorTime is a duration not a rate.
iorValue = 5;
saccScale = 1;
theta_saccStart = 0.2; %Saccade start threshold


% These parameters shouldn't change between models but do

% One thing to note about the field model is that field to discrete values
% should be calculated not by the magnitude of the location of the discrete
% value in the projecting field but by integrate the output field by the
% range of its kernel, weighted by the shape of that kernel.


if fieldModel
    sigma_ff = 1;
    sigma_fv = 1;
    sigma_cc = 1;
    sigma_prefpref =2;
    q_s = 0.5 * q;
    c_ff_7_83 = 5;
    c_ff_inh_7_81 = 2;
    c_cc_exc_8_99 =0.9;
    c_cc_inh_8_100 = 0.5;
    c_fv_7_82 = 20;
    c_gx_4_65 = 0.7;
    c_rx_6_76 = -5;
    c_ax_2_20 = 0.01;
    c_xr_5_70 = -0.5; %making this number larger makes the fixation force stronger.
    c_av_2_45 = 5;
    c_gg_4_66 = 0.2;
    c_prefgain_9_105 = 0.1;
    c_xcategory_changes = 10/(catFieldSize/30);
    c_xfeature_changes = 25/(featureFieldSize/60);
    advanceTrialThreshold = 0.95;
    CategoryChoiceThreshold = 0.97;
    c_clickimp_10_114 = 0.055;
    feedback_c_strength = 5;
    c_sa_3_53 = 8;
    c_as_2_41 = 1;
    c_aa_exc_2_34 = 0.05;
    c_va_1_5 = 0.0;
    sigma_aa_exc_2_35 = 10*gridRefinement;
    sigma_aa_inh_2_37 = 14*gridRefinement;
    c_vpref_1_13 = 0.5;
    c_vv_exc_1_6 = 6;
    c_va_gi_1_12 = 0.0;
    c_aa_inh_2_36 = 0;
    c_rg_6_77 = 6;
elseif ismember(structure,[3 7]) && ~fieldModel
    c_cf_8_98 = 1;
    c_clickc_10_113 = 2;
    c_sa_3_53 = 12;
    feedback_c_strength = 0.7;
    c_prefgain_9_105 = 0.2;
end

% This is used to supress foveal input to the saccade field from the
% attention field.
fovSuppression = c_sfov_3_56 * (1 - gauss2d_LAG1(visualHalfSize, sigma_sfov_3_55));

% Right now the assumption is that the fields are about 10 times the size of
% the neuron models. This should be updated at a later point.

if fieldModel
    kernelWidthMultiplier = 40*gridRefinement;
else
    kernelWidthMultiplier = 4*gridRefinement;
end



%% 6. Field, neuron and weight initializations

% initialize various fields and neurons to their resting states.
field_v = zeros(visualFieldSize,visualFieldSize,visualLayers) + h_v_1;
field_a = zeros(spatialFieldSize, spatialFieldSize) + h_a_18;
field_s = zeros(spatialFieldSize, spatialFieldSize) + h_s_47;
neuron_r = h_r_75;
neuron_x = h_x_68;
neuron_g = h_g_63;
%neuron_i = h_i;
neuron_click = h_click_107;
neuron_imp = h_imp_122;
input_v = zeros(visualFieldSize,visualFieldSize,visualLayers); %+1 for the clear added things
input_a = zeros(spatialFieldSize, spatialFieldSize);
activation_threshold_features = h_f_80+1;
activation_threshold_category = 0;
ior_input = zeros(spatialFieldSize,spatialFieldSize);
neuron_fbButton = h_f_80bButton;
neuron_fbButtonExp = h_f_80bButtonExp;


% In order to put things on the visual field there are three options. 1)
% Manually set colour components in experiment constructor. 2) Add input
% directly to the clear layer. 3) Break an image down using a
% spectral2slice function.


% wT := associative weight matrix. For example, this might be 4x6 because there are 6 feature
%neurons and 4 categories. One weight to each category from each feature.
if fieldModel
    field_f = zeros(featureFieldSize, 1) + h_f_80; %number of features.
    field_pref = zeros(featureFieldSize, 1) + h_pref_106; %number of features.
    field_c = zeros(catFieldSize, 1) + h_c_97; %number of categories.
    wT = 0.001*rand(catFieldSize,featureFieldSize) + h_wT + 0.9;
else
    neurons_f = zeros(featureNum, 1) + h_f_80; %number of features.
    neurons_pref = zeros(featureNum, 1) + h_pref_106; %number of features.
    neurons_c = zeros(catNum, 1) + h_c_97; %number of categories.
    wT = 0.001*rand(catNum,featureNum) + h_wT + 0.9;
end

output_wT_old = sigmoid(wT, beta_wT, inflec_wT);

% Setting the convolution kernels.

% a word about what's happening here. These tuples of code, specify 1. the
% size of the excitatory kernel input based on the smallest value of the
% arguments listed and 2. a normalized gaussian disitribution given a
% width, a mean, and a variance. These gaussian kernals will presumably be
% used in the convolution of the field they are included with. So, the vv
% kernel will be used in the convolution that determines the activations on
% a particular field.

% gaussNorm2d makes sure that width of kernel_vv_exc_spt is at least
% 2*kSize_vv_exc_spt +1.
kSize_vv_exc_spt = min(round(kernelWidthMultiplier * sigma_vv_exc_spt_1_7),visualFieldSize);
kernel_vv_exc_spt = gaussNorm2d(kSize_vv_exc_spt,sigma_vv_exc_spt_1_7);

kSize_vv_inh_s_47pt = min(round(kernelWidthMultiplier * sigma_vv_inh_s_47pt_1_9),visualFieldSize);
kernel_vv_inh_s_47pt = gaussNorm2d(kSize_vv_inh_s_47pt, sigma_vv_inh_s_47pt_1_9);

% Notice how fields that have self excitatory and self inhibitory
% connections are activated as a function of the difference between the
% two.

kSize_aa = min(round(kernelWidthMultiplier * max(sigma_aa_exc_2_35, sigma_aa_inh_2_37)), spatialFieldSize);
kernel_aa = c_aa_exc_2_34 * gaussNorm2d(kSize_aa, sigma_aa_exc_2_35) ...
    - c_aa_inh_2_36 * gaussNorm2d(kSize_aa, sigma_aa_inh_2_37);
kSize_ss = min(round(kernelWidthMultiplier * max(sigma_ss_exc_3_49, sigma_ss_inh_3_61)), spatialFieldSize);
kernel_ss = c_ss_exc_3_59 * gaussNorm2d(kSize_ss, sigma_ss_exc_3_49) ...
    - c_ss_inh_3_60 * gaussNorm2d(kSize_ss, sigma_ss_inh_3_61);

kSize_av = min(round(kernelWidthMultiplier * sigma_av_2_46), spatialFieldSize);
kernel_av = c_av_2_45 * gaussNorm2d(kSize_av, sigma_av_2_46);

kSize_va = min(round(kernelWidthMultiplier * sigma_va_1_4), spatialFieldSize);
kernel_va = c_va_1_5 * gaussNorm2d(kSize_va, sigma_va_1_4);



if fieldModel
    
    kSize_vpref = min(round(kernelWidthMultiplier * sigma_vpref), featureFieldSize);
    kernel_vpref = c_vpref_1_13 * gaussNorm(1:kSize_vpref,round(kSize_vpref/2), sigma_vpref);
    
    kSize_ff = min(round(kernelWidthMultiplier * sigma_ff), featureFieldSize);
    kernel_ff = c_ff_7_83 * gaussNorm(1:kSize_ff,round(kSize_ff/2), sigma_ff);
    
    kSize_fv = min(round(kernelWidthMultiplier * sigma_fv), featureFieldSize);
    kernel_fv = c_fv_7_82 * gaussNorm(1:kSize_fv,round(kSize_fv/2), sigma_fv);
    
    kSize_prefpref = min(round(kernelWidthMultiplier * sigma_prefpref), featureFieldSize);
    kernel_prefpref = c_prefpref_exc_9_110 * gaussNorm(1:kSize_prefpref,round(kSize_prefpref/2), sigma_prefpref);
    
    kSize_cc = min(round(kernelWidthMultiplier * sigma_cc), featureFieldSize);
    kernel_cc = c_cc_exc_8_99 * gaussNorm(1:kSize_cc,round(kSize_cc/2), sigma_cc);
    
else
    kSize_vpref = min(round(kernelWidthMultiplier * sigma_vpref), featureNum);
    kernel_vpref = c_vpref_1_13 * gaussNorm2d(kSize_vpref, sigma_vpref);
end


kSize_sa = min(round(kernelWidthMultiplier * sigma_sa_3_54), spatialFieldSize);
kernel_sa = c_sa_3_53 * gaussNorm2d(kSize_sa, sigma_sa_3_54);

kSize_as = min(round(kernelWidthMultiplier * sigma_as_2_42), spatialFieldSize);
kernel_as = c_as_2_41 * gaussNorm2d(kSize_as, sigma_as_2_42);

kSize_q_spatial = min(round(kernelWidthMultiplier * sigma_q_spatial_17), spatialFieldSize);
kernel_q_spatial = gaussNorm2d(kSize_q_spatial, sigma_q_spatial_17);



if length(varargin)>2
    
    for i=1:length(varargin{3})
        eval([varargin{3}{i} ' = ' num2str(varargin{4}(i))])
    end
    
    trialLogger(1,22) = {'Parameters'};
    trialLogger(2,22) = {varargin{3}};
    trialLogger(1,23) = {'Model_Number'};
    trialLogger(2,23) = {subjectNumber};
    
end


%% 8. Run-time parameters

breakOnFirstSaccade = false;  %for debugging.
tMax = 30000; %maximum possible time steps.

%% 9. Visualization
% GUI button assignment.
PLAY = 1;
PAUSE = 2;
STOP = 3;
QUIT = 4;

% Intializing GUI state.
state = PLAY;
initExperiment = true; % starts initialize experiment.
initTrial = true;
neuronHistory = [];

if visualize || startVisualize
    
    
    hFig = openfig('LAG1UI.fig');
    hAllAxes = findobj(gcf,'type','axes')
    
    Visual_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Visual Field'),get(hAllAxes,'Title')))).Children;
    Attention_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Spatial Attention Field'),get(hAllAxes,'Title')))).Children;
    Saccade_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Saccade Motor Field'),get(hAllAxes,'Title')))).Children;
    
    hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Visual Field'),get(hAllAxes,'Title')))).CLim=[0 1];
    hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Spatial Attention Field'),get(hAllAxes,'Title')))).CLim=[0 1];
    hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Saccade Motor Field'),get(hAllAxes,'Title')))).CLim=[0 1];
    
    Gaze_Change_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Gaze Change'),get(hAllAxes,'Title')))).Children;
    Fixation_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Fixation Excitation'),get(hAllAxes,'Title')))).Children;
    Saccade_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Saccade Reset'),get(hAllAxes,'Title')))).Children;
    FB_Exp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'FB  Expectation'),get(hAllAxes,'Title')))).Children;
    FB_Det_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'FB detection'),get(hAllAxes,'Title')))).Children;    
    Category_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Category Neurons'),get(hAllAxes,'Title')))).Children;
    Ft_Exp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Feature Expectations'),get(hAllAxes,'Title')))).Children;
    Ft_Det_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Feature Detection'),get(hAllAxes,'Title')))).Children;
    Motor_Click_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Motor Click Neuron'),get(hAllAxes,'Title')))).Children;
    Trial_Imp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Trial Impatience'),get(hAllAxes,'Title')))).Children;
    Trial_Num_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children(1);
    Exp_Time_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children(2);
    F1_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children(3);
    F2_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children(4);
    F3_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children(5);
    Max_Att_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Max Attention Activity'),get(hAllAxes,'Title'))));
    Weight_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Weight Matrix'),get(hAllAxes,'Title')))).Children;
    Max_Sacc_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Max Saccade Activity'),get(hAllAxes,'Title'))));
    Screen_Display_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Screen w/fovea'),get(hAllAxes,'Title')))).Children;
    Raw_Max_Att_Position = Max_Att_Axis.Position; % position of first axes
    Raw_Max_Att_Axis = axes('Position',Raw_Max_Att_Position,'XAxisLocation','top',...
        'YAxisLocation','left',...
        'Color','none','ylim',[-5.5 0],'YTick',[],'xlim',[1 51],'XTick',[]);
    
    Raw_Max_Sacc_Position = Max_Sacc_Axis.Position; % position of first axes
    Raw_Max_Sacc_Axis = axes('Position',Raw_Max_Sacc_Position,...
        'XAxisLocation','top',...
        'YAxisLocation','left',...
        'Color','none','ylim',[-6 13],'YTick',[],'xlim',[1 51],'XTick',[]);
    
    Gain_Applied_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
        'Gain'),get(hAllAxes,'Title')))).Children;
    
    Spatial_Gain_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,'Space Coded'),get(hAllAxes,'Title')))).Children;
    
    set(hAllAxes(find(cellfun(@(x) strcmp(x.String,'Feature Detection'),get(hAllAxes,'Title')))),'YLim',[1 featureNum]);
    set(hAllAxes(find(cellfun(@(x) strcmp(x.String,'Feature Expectations'),get(hAllAxes,'Title')))),'YLim',[1 featureNum]);
    set(hAllAxes(find(cellfun(@(x) strcmp(x.String,'Gain'),get(hAllAxes,'Title')))),'YLim',[1 featureNum]);
    set(hAllAxes(find(cellfun(@(x) strcmp(x.String,'Weight Matrix'),get(hAllAxes,'Title')))),'YLim',[1 featureNum],'XLim',[1 catNum]);
    set(hAllAxes(find(cellfun(@(x) strcmp(x.String,'Category Neurons'),get(hAllAxes,'Title')))),'YLim',[0.5 catNum+0.5]);
        
    %set(F3_circle, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 1 1], 'MarkerSize', 10);hold on
    
    
    
end

environLoaded =1; %Set to 0 to load a specific trial state.


while state ~= QUIT
    if initExperiment
        initExperiment = false;
        if profiling
            profile on
        end
    end
    
    
    
    
    while state == PLAY && trialNum <= maxTrials
        
        %% 10. TRIAL INITIALIZATION
        tic
        if initTrial
            
            
            % Neuron initializations
            field_v(:) = h_v_1;
            field_a(:) = h_a_18;
            field_s(:) = h_s_47;
            if fieldModel
                field_pref(:) = h_pref_106;
                field_f(:) = h_f_80;
                field_c(:) = h_c_97;
            else
                neurons_pref(:) = h_pref_106;
                neurons_f(:) = h_f_80;
                neurons_c(:) = h_c_97;
            end
            neuron_r(:) = h_r_75;
            neuron_x(:) = h_x_68;
            neuron_g(:) = h_g_63;
            %neuron_i(:) = h_i;
            ior_input = zeros(spatialFieldSize,spatialFieldSize);
            input_v(:) = 0;
            input_a(:) = 0;
            neuron_click = h_click_107;
            neuron_imp = h_imp_122;
            neuron_fbButton = h_f_80bButton;
            neuron_fbButtonExp = h_f_80bButtonExp;
            
            c_cf = c_cf_8_98; %In order to prevent category activation during feedback
            firstFBTrialFlag=1; % 1 if this is the first feedback timestep in phase 4 (feedback viewing). Gets turned off on the 2nd p4 fb timestep.
            StimulusOnset =0; % This records the time to first fixation in TET time units. [Blair lab standard]
            trialPhase = 1; % This specifies the fixation cross component of the trial. [Blair lab standard]
            % trialPhase means: 1 = fixation cross; 2 = stimulus presentation; 3 = ; 4 = feedback presentation
            FixationOnset = TETTime; %This specifies the absolute start of the trial.
            kernel_av = c_av_2_45 * gaussNorm2d(kSize_av, sigma_av_2_46); %In case of ending on saccade
            decisionFlag = 0; % have we chosen category yet in this trial
            fbStartFlag = 0; %Unlike firstFBTrialFlag, this indicates that generally we're in the feedback phase and stays on.
            fbStartTime = 0; %Starts counting as soon as the click neuron reaches threshold - not when feedback is actually present.
            fixationImp = 0; %Records the time step that the extant fixation started at.
            spatialFieldFovea = [round(spatialFieldSize/2) round(spatialFieldSize/2)]; %The center of the field in matrix coordinates.
            saccadeTimes = []; %variable length vector of the saccade durations for a trial
            outstandingGainWeight = 0; %a scalar estimating how much important stuff remains to be seen
            feedbackPhase = 0; %set to 1 precisely when phase 4 begins, perhaps redundant.
            output_c = 0; % in order to calculate a change value between a category change input to the fixation neurons we need this initailized
            output_f = 0; % in order to calculate a change value between a feature change input to the fixation neurons we need this initailized
            output_fbButton = 0; % in order to calculate a change value between a fb feature change input to the fixation neurons we need this initailized
            trialFixations = []; % Records the features fixated over the course of the trial.
            trialFixationsFB = []; %prints feedback fixation locations as a debugging
            trialFix = nan;
            fixationChanges = [];
            fixationChangesRecord = [];
            ColourOnset = 0;
            valueGain = 0*wT;
            featureValue = 0;
            fbViewed = 0;
            featureDetected = 0;
            feedback_c = zeros(catFieldSize,1);
            fbButtonDetection = 0;
            c_vFBnovelty_1_4 = 10;
            Response = 0; %Detect if a response has yet been made.
            
            if (transfer.switch == 1) || (transfer.switch == 0) %i.e. training phase.
                resultsStimulus = [resultsStimulus;{stimulus(trialNum, : )}]; %Records the stimulus values.
                CorrectResponse = [CorrectResponse;{category(trialNum)}]; % Records the associated category value.
                stimFeature_v = stimulus(trialNum,:); %Sort of redundant but these ARE the stimulus values now, not just the record of them.
                binaryFeatureValues = stimFeature_v;
                sceneStimColors = colours(trialNum, :); %Get the colour values for the trial.
            elseif transfer.switch == 2 %i.e. transfer phase
                resultsStimulus = [resultsStimulus;{stimulus(trialNum-1000, : )}]; %Records the stimulus values.
                CorrectResponse = [CorrectResponse;{category(trialNum-1000)}]; % Records the associated category value.
                stimFeature_v = stimulus(trialNum-1000,:); %Sort of redundant but these ARE the stimulus values now, not just the record of them.
                binaryFeatureValues = stimFeature_v;
                sceneStimColors = colours(trialNum-1000, :); %Get the colour values for the trial.
            end
            
            
            % This allows counterbalanced features to get appropriately
            % described in terms of their feature value. Note: this is now
            % done much more simply in counterBalanceChecker.m
            
            binaryFeatureValues(binaryFeatureValues==1)=0;
            binaryFeatureValues(binaryFeatureValues==2)=1;
            binaryFeatureValues(binaryFeatureValues==3)=0;
            binaryFeatureValues(binaryFeatureValues==4)=1;
            binaryFeatureValues(binaryFeatureValues==5)=0;
            binaryFeatureValues(binaryFeatureValues==6)=1;
            
            if featureNum==8
                binaryFeatureValues(binaryFeatureValues==7)=0;
                binaryFeatureValues(binaryFeatureValues==8)=1;
            end
            %Begins visualization on a pre-specified trial. This is a
            %debugging helper.
            if trialNum == startVisualize
                visualize = 1;
                beep;pause(2);beep;pause(1.5);beep;pause(1.1);beep;pause(.8);beep;pause(.7)
            end
            
            if trialNum == endVisualize
                visualize = 0;
            end
            
            
            fixCrossInd = randi(3); %This is currently unused but allows for a random fixation cross location.
            gazeShift=[0 0]*gridRefinement; % Where is gaze in relation to the spatial origin. Needed for movement calculation.
            slidingVisX = 1+visualHalfSize+(-spatialHalfSize:spatialHalfSize)-gazeShift(1); %Re-orients the spatial field in relation to the visual field on every eye movement.
            slidingVisY = 1+visualHalfSize+(-spatialHalfSize:spatialHalfSize)-gazeShift(2); %ditto but for the y axis.
            stimSize_v = targetSize*ones(1,locationNum); % Specifies how large each feature should be.
            stimTimes_v = repmat([1, tMax],totalLocations,1); % How long should each stim be presented for on a particular trial. i.e. 30000 timesteps.
            xshift = 0; %Simply initalizing
            yshift = 0; %Simply initalizing
            input_i = zeros(tMax,1);  %Simply initalizing
            stimuli_v = zeros(spatialFieldSize, spatialFieldSize, featureNum); %Initializing the stimulus presentation structure
            stimWeights_v = zeros(tMax, totalLocations); % Specifies a salience bias to a particular feature.
            sceneStimActive = zeros(tMax, totalLocations); %Initializes a large matrix specifying whether a feature is active at a particular timestep
            saccStartTimes = [];
            saccEndTimes = [];
            saccadeInProgress = false;
            gazeChangeCounter = 0;
            
            
            
            %% Making the stimuli appear on the visual field
            for i = 1 : totalLocations
                
                
                if i < totalLocations
                    
                    %Input to the visual field.
                    [YY,XX]=meshgrid(-spatialHalfSize:spatialHalfSize);
                    % important point: We want the first coordinate to be x, and
                    % the second y. When we plot these fields later we will have to
                    % plot transpose to make it look right, because matlab plots
                    % first coordinate as vertical. This is why YY and XX are
                    % reversed to what you'd expect in the above line.
                    RR=sqrt( (stimPos_x(i)-XX).^2 + (stimPos_y(i)-YY).^2); % distances from feature
                    spatialInput=(RR<stimSize_v(i)/2); %Puts a 1 on the field wherever a feature is
                    
                    if stimFeature_v(i) ~= -1
                        stimuli_v(:, :, stimFeature_v(i)) = spatialInput; %Associates the feature existence matrix with the right feature value in the 3rd dimension of the visual field.     
                    end
                    stimWeights_v(stimTimes_v(i, 1):stimTimes_v(i, 2), i) = 1;
                    featureColours = []; 
                    %This will create a 3 x n matrix with n being the number of features and 3 being the components of R-G-B
                    
                    for k = 1:locationNum %3 being the components of RGB. 3 rows X featureNum columns

                        if sceneStimColors(k)==-1
                            featureColours = [featureColours [0;0;0]];
                        else
                            featureColours = [featureColours [[sceneColors(1, sceneStimColors(k), 1)];[sceneColors(1, sceneStimColors(k), 2)];[sceneColors(1, sceneStimColors(k), 3)]]];
                            
                        end

                    end
                    
                    %Keeps these stims active for the prespecified max trial time.
                    sceneStimActive(stimTimes_v(i, 1):stimTimes_v(i, 2), i) = 1;
                else
                    [YY,XX]=meshgrid(-spatialHalfSize:spatialHalfSize);
                    RR=sqrt((stimPos_x(i)-XX).^2 + (stimPos_y(i)-YY).^2);
                    spatialInput=(RR<fbButtonSize/2);
                    stimuli_v(:, :, visualLayers) = spatialInput*c_inp_1;
                    stimWeights_v(stimTimes_v(i, 1):stimTimes_v(i, 2), i) = 1;
                    sceneStimActive(stimTimes_v(i, 1):stimTimes_v(i, 2), i) = 1;
                    
                end
                
            end
            
            %%
            
            % clear visual field and fill the right portion of it with
            % the simulus.
            shiftedStimuli_v=zeros(visualFieldSize,visualFieldSize,visualLayers);
            shiftedStimuli_v(slidingVisX,slidingVisY,:)=stimuli_v;
            
            
            trialTimeStep = 1;
            initTrial = false;
            display('##############################################')
            display(['This is trial number ' num2str(trialNum) ' of ' num2str(maxTrials) ' total.'])
        end
        
        
        %% 11. START TRIAL
        while state == PLAY &&  trialTimeStep <= tMax
            
            
            
            %Output from the different fields is determined after sigmoiding
            %the current activations of those fields.
            output_cOld = output_c;
            output_fOld = output_f;
            output_fbButtonOld = output_fbButton;
            output_v = sigmoid(field_v, beta_v_10, inflec_v_11);
            output_a = sigmoid(field_a, beta_a_32, inflec_a_33);
            output_s = sigmoid(field_s, beta_s_43, inflec_s_44);
            if fieldModel
                output_f = sigmoid(field_f, beta_f_85, inflec_f_86);
                output_pref = sigmoid(field_pref, beta_pref_15, inflec_pref_16);
                output_c = sigmoid(field_c, beta_c_101, inflec_c_102);
            else
                output_f = sigmoid(neurons_f, beta_f_85, inflec_f_86);
                output_pref = sigmoid(neurons_pref, beta_pref_15, inflec_pref_16);
                output_c = sigmoid(neurons_c, beta_c_101, inflec_c_102);
            end
            output_r = sigmoid(neuron_r, beta_r_26, inflec_r_27);
            output_x = sigmoid(neuron_x, beta_x_50, inflec_x_51);
            output_g = sigmoid(neuron_g, beta_g_30, inflec_g_31);
            %output_i = sigmoid(neuron_i, beta_i, inflec_i);
            output_a_fovSup = output_a .* fovSuppression(slidingVisX,slidingVisY);
            output_wT = sigmoid(wT, beta_wT, inflec_wT);
            output_click = sigmoid(neuron_click, beta_click_116, inflec_click_117);
            output_imp = sigmoid(neuron_imp, beta_imp_120, inflec_imp_121);
            output_fbButton = sigmoid(neuron_fbButton, beta_f_85bButton, inflec_f_86bButton);
            output_fbButtonExp = sigmoid(neuron_fbButtonExp, beta_f_85bButtonExp_23, inflec_f_86bButtonExp_24);
            
            
            %The trial decision threshold for phase 2.
            %decisionFlag makes it so we're only here once.
            
            %single click
            if (output_click > CategoryChoiceThreshold) && (output_click < advanceTrialThreshold)  && (decisionFlag == 0) && (trialPhase == 2) && ~saccadeInProgress
                decisionFlag = 1; %Cannot come back here on this trial.
                fbStartFlag = 1; %The person has clicked but they have not looked away, they are still in phase 2 as far as we're concerned.
                output_wT_old = output_wT;
                ColourOnset = TETTime; % This is a phase 3 thing, unneeded.
                feedbackPhase = 1;
                
                %double click or end of feedback
                %This is the feedback threshold that ends a trial.
            elseif (output_click > advanceTrialThreshold) && trialPhase==4 && ~saccadeInProgress
                decisionFlag = 1;
                trialLogger{2,8} = [trialLogger{2,8};trialTimeStep - fbStartTime];
                trialLogger{2,9} = [trialLogger{2,9};{output_wT - output_wT_old}];
                saccStartTimes = [saccStartTimes; trialTimeStep];
                saccEndTimes = [saccEndTimes; trialTimeStep];
                saccadeTimes = [saccadeTimes 0];
                trialTimeStep = tMax;
                
                %LAG1 should only come here when there's been a
                %doubleclick.
                if ~Response
                    Response = randsample(catNum,1,'true',softermax(output_c',temperature));
                    FeedbackRT = 0;
                    
                    %Is there a tie? Randomly select a winner.
                    if length(Response)>1
                        Response = Response(randi(length(Response)));
                    end
                    %Correct?
                    if Response == CorrectResponse{end}
                        TrialAccuracy=1;
                    else
                        TrialAccuracy=0;
                    end
                end
                
                
            end
            
            %This catches bad behaviour and ends the simulation.
            if trialTimeStep==badTimestepNum || length(trialFixations) >= badFixationNum
                
                if optimizationModeFlag
                    endFit=NaN;
                    return
                end
                trialLogger{2,8} = [trialLogger{2,8};trialTimeStep - fbStartTime];
                trialLogger{2,9} = [trialLogger{2,9};{[]}];
                trialTimeStep = tMax;
                
            end
            
            %% 12. Feedback phase
            if decisionFlag && ~feedbackOff
                if firstFBTrialFlag == 1 % If this is the first timestep of the feedback

                    %Greedy or softMax?
                    if softMax
                        Response = randsample(catNum,1,'true',softermax(output_c(catIndex(:))',temperature));
                        %Loop through temperature values to get differing
                        %accuracy curves
                        if ~antiHebbCorr
                            accVec =[];
                            for SoftmaxTemperature = accuracyLevels{1}'
                                accVec = [accVec randsample(catNum,1,'true',softermax(output_c',SoftmaxTemperature))==CorrectResponse{end}];
                            end
                            accuracyLevels{1,1} = [accuracyLevels{1,1} [accVec]'];
                        end
                    else
                        Response = find(max(output_c) == output_c)
                    end
                    
                    fbStartTime = trialTimeStep; %Another feedback timer. This time starting on the first time step of phase 4.
                    
                    %Is there a tie? Randomly select a winner.
                    if length(Response)>1
                        Response = Response(randi(length(Response)));
                    end
                    
                    %Correct?
                    if Response == CorrectResponse{end}
                        TrialAccuracy=1;
                    else
                        TrialAccuracy=0;
                    end
                    trialLogger{2,14} = [trialLogger{2,14}; CorrectResponse{end}];
                    trialLogger{2,15} = [trialLogger{2,15}; Response];
                    trialLogger{2,18} = [trialLogger{2,18}; TrialAccuracy];
                    %Log the category and feature activations
                    trialLogger(2,1:4) = [{[trialLogger{2,1} neurons_c]} {[trialLogger{2,2} output_c]} {[trialLogger{2,3} neurons_f]} {[trialLogger{2,4} output_f]}];
                    trialLogger(2,5:6) = [{[trialLogger{2,15}; {wT}]} {[trialLogger{2,16}; {output_wT}]}];
                    
                    firstFBTrialFlag = 0;
                end
                
                %% 13. Hebbian learning.
                
                %Checks that there is coactivation of feature and category
                %neurons.
                
                if fbViewed
                    if fieldModel && ~firstFBTrialFlag
                        
                        feedback_c = 0; %initialize the feedback input.
                        if TrialAccuracy || (feedbackType~=3) && (feedbackType~=4)
                            % This gets applied in the field update.
                            feedback_c = output_fbButton*feedback_c_strength*gauss(1:catFieldSize,catIndex(CorrectResponse{end}),sigma_cc);
                        else %come in here for unsupervised learning.
                            field_c = output_fbButton*feedback_c_strength*gauss(1:catFieldSize,catIndex([1:Response-1 Response+1:end]),sigma_cc);
                        end
                        
                    elseif ~fieldModel && ~firstFBTrialFlag
                        
                        if TrialAccuracy || (feedbackType~=3) && (feedbackType~=4)
                            feedback_c = -output_fbButton*feedback_c_strength*ones(catNum,1);
                            feedback_c(CorrectResponse{end}) = output_fbButton*feedback_c_strength;
                        else %come in here for unsupervised learning.
                            neurons_c([1:Response-1 Response+1:end]) = output_fbButton*feedback_c_strength;
                        end
                        
                    end
                    
                    firstFBTrialFlag = 0; % don't come back here on this trial.
                    
                    if fieldModel
                        wT = associator(wT, tau, field_f, field_c, ...
                            learning_rate, TrialAccuracy, ...
                            feedbackType,activation_threshold_features, ...
                            activation_threshold_category,antiHebbCorr, ...
                            decay_rate,deltaT,antiHebbRate,Response);
                    else
                        
                        %% Associative learning
                        % Document reference: \ref{eq:hebbEqs}
                        % $$\tau \dot w(i,j,t) = -w(i,j,t) + {h_1}$$
                        %
                        
                        wT = associator(wT, tau, deltaT, neurons_f, neurons_c, ...
                            learning_rate, TrialAccuracy, feedbackType,...
                            activation_threshold_features,activation_threshold_category, ...
                            antiHebbCorr,decay_rate,antiHebbRate,Response);
                        
                    end
                end
                
            end
            %sigmoids the integrated foveal feature values.
            
            %The presence of information of a particular value
            %currently being foveated is represented simply as on or off,
            %in the size it takes up within the fovea. We sum this to get a
            %value to feed into the feature neurons.
            
            %featureDetection_range=(visualHalfSize-fovealHalfSize):(visualHalfSize+fovealHalfSize); % range of indices for feature detection
            
            fovealMask = repmat(gaussNorm2d(visualHalfSize, sigma_fov_mask),[1,1,featureNum]);
            featureDetection_exc = sigmoid(sum(sum(output_v(:,:,1:featureNum).*fovealMask)),beta_ft_input_87,inflec_ft_input_88);
            
            if fbViewed
                fovealMask = repmat(gaussNorm2d(visualHalfSize, sigma_fov_mask),[1,1,visualLayers]);
                fbButtonDetection = sigmoid(sum(sum(output_v(:,:,1:visualLayers).*fovealMask)),beta_ft_input_87,inflec_ft_input_88);
                fbButtonDetection = fbButtonDetection(:,:,end); %This is was accidentally set to 7 but makes no difference. It should be fixed in version 2.0
            end
            
            if fieldModel
                % "Detect colour" In the future this will have to be thought
                % out more. The choice of colour representation doesn't matter
                % at all as it stands. i.e. whether the colour fades at the
                % spatial edges or is enhanced by contrast or turns into a
                % Gaussian between the visual field and the feature field.
                % Right now it's the latter.
                input_fv = 0;
                for i = 1:featureNum
                    
                    %When the pairings aren't known, this will be needed:
                    
                    input_fv = input_fv + (featureDetection_exc(i)*gauss([1:featureFieldSize],featureAntiCorrelations(i,2),2) - featureDetection_exc(featureAntiCorrelations(i,3))*gauss([1:featureFieldSize],featureAntiCorrelations(i,2),2));
                    
                end
                input_fv = c_fv_7_82 * input_fv;
                
            else
                %Squeeze removes the singleton dimensions which exist
                %because we're only interested in the values of the 3rd
                %dimension of the visual field.
                featureDetection_exc = squeeze(featureDetection_exc);
                %featureAntiCorrelations(:,3) records the anti-correlated
                %feature values of the experiment. This will be more useful
                %if we're using a feature field instead of feature neurons.
                featureDetection_inh = -featureDetection_exc(featureAntiCorrelations(:,3)');

                %input from the visual field to feature neurons is a function
                %of what has been detected and what logically cannot be,
                %based on the instructions.
                input_fv = c_fv_7_82 * (featureDetection_exc + featureDetection_inh);
            end
            
            %% 14. Saccade and rotation mechanics.
            
            % determine start of saccade and suppress input. output_r must
            % cross a threshold for the saccade to begin.
            if (~saccadeInProgress && output_r >= theta_saccStart && trialTimeStep~=tMax)
                
                fixationImp = trialTimeStep; %fixation impatience timer reference
                saccStart = trialTimeStep;
                fixationChangesRecord = [fixationChangesRecord;sum(fixationChanges)];
                fixationChanges = [];
                
                
                if decisionFlag
                    trialPhase=4;
                    StimulusRT = TETTime - StimulusOnset; % phase 2 timing.
                end
                
                
                
                featureDetected = 0;
                
                %Increment gaze change counter
                gazeChangeCounter = gazeChangeCounter + 1;
                %Record the start time and cardinality of a saccade.
                saccStartTimes(gazeChangeCounter, 1) = trialTimeStep; %record the time a saccade begins.
                
                saccadeInProgress = true;
                kernel_av = 0 * gaussNorm2d(kSize_av,sigma_av_2_46); %kill input from the visual field into the attention field.
                
                
                [max_output_s,ind]=max(output_s(:)); %find the the largest peak on the saccade field.
                [target_gaze_shift1,target_gaze_shift2] =ind2sub(size(output_s),ind); % find the location of of that maximum
                ior1 = spatialFieldFovea(1); %Set the IOR location to the extant coordinates of the fovea.
                ior2 = spatialFieldFovea(2); %ditto for Y
                xdiff = target_gaze_shift1-spatialFieldFovea(1); %Get the magnitude of the expected X change
                ydiff = target_gaze_shift2-spatialFieldFovea(2); %ditto for y
                saccadeDistance = hypot(xdiff, ydiff);
                %larger saccScale = shorter saccades.
                saccadeDuration = (saccadeDistance/saccScale);
                
                
                xshift = xdiff/saccadeDuration; %units to travel/timestep X
                yshift = ydiff/saccadeDuration; %ditto for Y.
                xError=0; %reset the error on the actual vs. estimated end point of the last saccade on X.
                yError=0; %reset the error on the actual vs. estimated end point of the last saccade on Y.
                
                
                
                %This will create an IOR to the fixated location. The stim
                %positions were changed to references, not absolute
                %values. This code will only produce an IOR on a feature and
                %not a location. This is something to consider for the future.
                %Currently the IOR is set to last for 3 seconds.
                createIOR = 0;
                for i=1:totalLocations
                    if ( norm(gazeShift-[ stimPos_x(i) stimPos_y(i)])<fovealHalfSize ) %only create IOR if something was fixated.
                        createIOR = 1;
                    end
                end
                if createIOR
                    input_i(trialTimeStep:trialTimeStep+iorTime) = iorValue;
                else
                    input_i(trialTimeStep:end) = 0;
                end
                c_ax1 = c_ax_2_20;
                c_ax_2_20=0;
                category_change_coeff_temp = c_xcategory_changes;
                c_xcategory_changes = 0;
                feature_change_coeff_temp = c_xfeature_changes;
                c_xfeature_changes = 0;
                fbButton_changes_temp = c_xfbButton_changes;
                c_xfbButton_changes = 0;
                
                
                if (output_click > advanceTrialThreshold)
                    decisionFlag = 1;
                    trialLogger{2,8} = [trialLogger{2,8};trialTimeStep - fbStartTime];
                    trialLogger{2,9} = [trialLogger{2,9};{output_wT - output_wT_old}];
                    saccEndTimes = [saccEndTimes; trialTimeStep+1];
                    saccadeTimes = [saccadeTimes 0];
                    trialTimeStep = tMax;
                    
                    if ~Response
                        Response = randsample(catNum,1,'true',softermax(output_c',temperature));
                        FeedbackRT = 0;
                        
                        %Is there a tie? Randomly select a winner.
                        if length(Response)>1
                            Response = Response(randi(length(Response)));
                        end
                        
                        %Correct?
                        if Response == CorrectResponse{end}
                            TrialAccuracy=1;
                        else
                            TrialAccuracy=0;
                        end
                    end
                    StimulusRT = TETTime - StimulusOnset; % phase 2 timing.
                end
                
                
                
            end
            
            
            %% Saccade travel
            
            if saccadeInProgress
                
                
                %Continually deal with the quantization error of the field
                xshift1 = xError + xshift;
                xError = xshift1 - floor(xshift1);
                yshift1 = yError + yshift;
                yError = yshift1 - floor(yshift1);
                
                %Update the current point of foveation on the spatial and
                %visual fields
                spatialFieldFovea(1) = spatialFieldFovea(1) + floor(xshift1);
                spatialFieldFovea(2) = spatialFieldFovea(2) + floor(yshift1);
                
                slidingVisX = slidingVisX-floor(xshift1);
                slidingVisY = slidingVisY-floor(yshift1);
                
                
                %Gazeshift is what actually moves the model.
                gazeShift(1) = gazeShift(1) +floor(xshift1);
                gazeShift(2) = gazeShift(2) +floor(yshift1);
                
                
                %% This should never be entered. There has been a saccade off screen.
                
                if (slidingVisX(1) < 1) || (slidingVisX(spatialFieldSize) > visualFieldSize) || (slidingVisY(1) < 1) || (slidingVisY(spatialFieldSize) > visualFieldSize)
                    subjectNumber=NaN; endFit=NaN;
                    error('There has been a saccade off screen.')
                end
                
                %%
                
                
                %Mapping visual stimuli into a retinal window. This
                %facilates interactions between attention and visual world.
                %Without this frame the input to v will be misaligned.
                shiftedStimuli_v = zeros(visualFieldSize,visualFieldSize,visualLayers);
                shiftedStimuli_v(slidingVisX,slidingVisY,:) = stimuli_v;
                
            end
            %%
            
            
            % Determine end of saccade and update input
            % Second feedbackType is needed to figure out when the saccade
            % actually ends.
            if (saccadeInProgress && (saccadeDuration + saccStart <= trialTimeStep)) && (trialTimeStep ~= tMax)
                
                
                saccadeTimes = [saccadeTimes trialTimeStep-saccStart;]; %Record the saccade time.
                saccEndTimes(gazeChangeCounter, 1) = trialTimeStep; %Record the saccade end time.
                click_ready=0; %Allows the model to enter phase 3 if needed.
                
                if (output_click > advanceTrialThreshold) && trialPhase==4
                    decisionFlag = 1;
                    trialLogger{2,8} = [trialLogger{2,8};trialTimeStep - fbStartTime];
                    trialLogger{2,9} = [trialLogger{2,9};{output_wT - output_wT_old}];
                    saccStartTimes = [saccStartTimes; trialTimeStep];
                    saccEndTimes = [saccEndTimes; trialTimeStep];
                    saccadeTimes = [saccadeTimes 0];
                    trialTimeStep = tMax;
                end
                
                if feedbackPhase && feedbackOff %This sets up the end of the trial if the model is doing transfer stims.
                    trialTimeStep=tMax;
                    Response = randsample(catNum,1,'true',softermax(output_c',temperature));
                    FeedbackRT = 0;
                    
                    %Is there a tie? Randomly select a winner.
                    if length(Response)>1
                        Response = Response(randi(length(Response)));
                    end
                    
                    TrialAccuracy=0;
                    
                end
                
                if trialPhase ==1
                    if StimulusOnset==0
                        StimulusOnset = TETTime;
                    end
                    trialPhase = 2;
                end
                
                if breakOnFirstSaccade
                    break;
                end
                
                %Create an IOR input.
                ior_input = (gauss2d_LAG1([spatialHalfSize spatialHalfSize],sigma_aIn_2_22,spatialHalfSize,spatialHalfSize,ior1-spatialHalfSize,ior2-spatialHalfSize)');
                
                
                %redfined this back to proper value after setting
                %kernel to 0 prior.
                kernel_av = c_av_2_45 * gaussNorm2d(kSize_av, sigma_av_2_46);
                
                % This set of conditions simply logs the feature viewed.
                
                for featureLocal=1:totalLocations
                    if ( norm(gazeShift-[ stimPos_x(featureLocal) stimPos_y(featureLocal)])<fovealHalfSize*2 )
                        trialFix=featureLocal; %Record the feature fixated.
                        if featureLocal == totalLocations %because the feedback button is the last one...
                            fbViewed = 1;
                            c_vFBnovelty_1_4 = 1; %Reduce the saliency of the fb button now.
                            c_aFBnovelty = 1;
                        end
                        if trialPhase ==2
                            trialFixations = [trialFixations,trialFix]; %Note which feature was looked at.
                        elseif trialPhase ==4
                            trialFixationsFB = [trialFixationsFB, trialFix];
                        end
                        featureDetected = 1;
                    elseif ~featureDetected && (featureLocal == totalLocations)
                        if trialPhase ==2
                            trialFixations = [trialFixations,totalLocations]; %Note the fb button was looked at.
                        elseif trialPhase ==4
                            trialFixationsFB = [trialFixationsFB, totalLocations];
                        end
                    end
                end
                
                saccadeInProgress = false;
                featureDetected = 0;
                
                %Apply the category and feature change input to the fixation
                %neuron if a feature is being foveated.
                c_xcategory_changes = category_change_coeff_temp;
                c_xfeature_changes = feature_change_coeff_temp;
                c_xfbButton_changes = fbButton_changes_temp;
                c_ax_2_20 = c_ax1;
                
            end
            
            
            %input_xChanges = c_xcategory_changes*sum(abs(output_c -output_cOld)),beta_ft_x_73,inflec_ft_x_74) + sigmoid(c_xfeature_changes*sum(abs(output_f -output_fOld)),beta_ft_x_73,inflec_ft_x_74);
            fixationChanges = [fixationChanges;c_xcategory_changes*sum(abs(output_c -output_cOld)) c_xfeature_changes*sum(abs(output_f -output_fOld)) c_xfbButton_changes*sum(abs(output_fbButton -output_fbButtonOld))];
            input_xChanges = sum(fixationChanges(end,:));
            
            
            %% 15. Neural update calculations.
            
            % Input to the visual field is modulated by whether there is a saccade happening or not. If there is, there is no input to the field.
            % The second part of this equation is what compacts the three
            % dimensions that were stored seperately to try and account for
            % different stimulus weights, down to one i.e. sum(var,3).
            
            % This c_vr_2 will suppress the visual field during a saccade.
            input_v = (1 - c_vr_1_3 * output_r) * c_vinp_1_2*shiftedStimuli_v;
            input_v(:,:,visualLayers) = c_vFBnovelty_1_4 * output_fbButtonExp * fbStartFlag * input_v(:,:,visualLayers);
            
            input_ar = c_ar_2_29*output_r*output_a; %Proportional to strength on attention field? weird claim. I dunno about this.
            
            output_f = output_f(:); %Just ensure it's a column vector.
            
            
            % Gained category activations
            if fieldModel

                featureValue(featureNum/2) = 0;
                for i = 1:2:featureNum
                    featureValue(i) = [sum(repmat(abs(output_wT(:,featureAntiCorrelations(i,1))-output_wT(:,featureAntiCorrelations(i+1,1))),1,2))];
                end
                
                valueGain = c_prefgain_9_105 * featureValue;
                valueGain = kron(valueGain,[ones(featureFieldSize/(featureNum))]);
                valueGain = repmat(valueGain,[catNum 1]);
                input_prefs = (output_wT' + valueGain') * output_c;
                input_prefc = c_prefc_9_109 * input_prefs;
                input_apref = c_apref_2_19 * output_pref;
                
            else
                
                %% Feature to Category gain
                % Document reference: \ref{eq:gainFunctionEq}
                % $$g(i,j,t) = w(i,j,t) + {c_{105}}\left| {w(i,j,t){\rm{  -  }}w(i,j',t)}\right|$$
                %
                
                
                if size(featureAntiCorrelations) > 0 % isolated locations for the feature values
                    
                    if ismember(structure, [1 2 4 5])
                        diffMatrix = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
                        valuesPerLocation = 2;
                        kronMultiple = ones(1,valuesPerLocation);
                    elseif ismember(structure, [3])
                        diffMatrix = [1 0 0 0;-1 0 0 0;0 1 0 0;0 -1 0 0;0 0 1 0;0 0 -1 0; 0 0 0 1;0 0 0 -1];
                        valuesPerLocation = 2;
                        kronMultiple = ones(1,valuesPerLocation);
                    elseif structure ==6 % Not in the blocking structure.
                        diffMatrix = [1 0 0 0 0;-1 0 0 0 0;0 1 0 0 0;0 -1 0 0 0;0 0 1 0 0;0 0 -1 0 0; 0 0 0 1 0;0 0 0 -1 0;0 0 0 0 1;0 0 0 0 -1];
                        valuesPerLocation = 2;
                        kronMultiple = ones(1,valuesPerLocation);
                    elseif ismember(structure, [7])
                        diffMatrix = [1 0 0 0;0 0 -1 0;0 1 0 0;0 0 0 -1];
                        valuesPerLocation = 1;
                        kronMultiple = ones(1,valuesPerLocation);
                    end
                    
                    featureValue = kron(output_c'*abs(output_wT*diffMatrix),kronMultiple);
                    
                    input_cf = c_cf * (c_cgain_8_105*(kron(abs(output_wT*diffMatrix),kronMultiple)) + output_wT) * (output_f);
                    
                    input_prefs = output_wT' * output_c;
                    input_prefs = input_prefs + c_prefgain_9_105*featureValue';  % why were there two 105s?
                    
                    
                else % feature values are not tied to location.
                    input_cf = c_cf * (output_wT) * (output_f); 
                    input_prefs = output_wT' * output_c; %NO GAIN.

                    
                end
                
            end
            
            input_prefc = c_prefc_9_109 * input_prefs; % category inputs to the preferences field
            input_apref = c_apref_2_19 * output_pref; % high-level attention to the attention field from the preference field.
            
            % As the most honest way to ascertain where the 3 features are
            % on the attention field is to use the inputs from the visual
            % field for the coordinates, this applies the preferences from
            % the individual feature values to each layer of the input from
            % the visual field as input directly to the attention field.
            
            attPrefTotal = 0;
            for i=1:featureNum
                attPrefTotal=attPrefTotal+ output_v(slidingVisX,slidingVisY,i)*input_apref(i);
            end
            
            attPref = conv2(attPrefTotal,kernel_av, 'same');
            
            
            %We can see activation of the feature hypothesis neurons get
            %converted to a field in the next few lines. f
            if fieldModel
                w_vpref = c_vpref_1_13*reshape(output_pref(featureAntiCorrelations(:,2)),1,1,featureNum);
            else
                w_vpref = c_vpref_1_13*reshape(output_pref,1,1,featureFieldSize);
            end
            
            
            
            %I'm not sure if the sum(output_pref) keeps the inhibition
            %proportional to the excitation or something?
            input_vpref = repmat(w_vpref,[visualFieldSize,visualFieldSize, 1]) - c_vpref_gi_1_14 * sum(output_pref);
            input_vpref(:,:,visualLayers) = 0;
            
            %The activation of each each dimension of the feature is summed
            %together. This creates essentially a preference for a particular
            %location, not a particular feature value.
            
            input_vv_exc=zeros(visualFieldSize,visualFieldSize,visualLayers);
            input_vv_inh=zeros(visualFieldSize,visualFieldSize,visualLayers);
            
            % This was massively confusing to come back to. It makes sense
            % after a while of thinking though, in that the inhibitory
            % component needs to come from all the space
            for i=1:visualLayers
                input_vv_exc(:,:,i) = conv2(output_v(:,:,i), kernel_vv_exc_spt, 'same');
                input_vv_inh(:,:,i) = conv2(sum(output_v,3), kernel_vv_inh_s_47pt, 'same');
            end
            input_vv = c_vv_exc_1_6 * input_vv_exc - c_vv_inh_1_8 * input_vv_inh;
            
            input_aFBbutton = c_afbButton_2_25 * c_aFBnovelty* output_fbButtonExp * output_a(spatialFBposition(1)+spatialHalfSize,spatialFBposition(2)+spatialHalfSize)' * gauss2d_LAG1([spatialHalfSize spatialHalfSize],sigma_av_2_46,spatialHalfSize,spatialHalfSize,spatialFBposition(1),spatialFBposition(2))';

            
            
            input_aa = conv2( output_a,kernel_aa, 'same') - c_aa_gi_2_38 * sum(output_a(:));% - input_aFBbutton_inh';
            input_ss = conv2( output_s,kernel_ss, 'same') - c_ss_gi_3_3_62 * sum(output_s(:));
            
            
            input_rr = c_rr_6_79 * output_r;
            input_av = conv2( sum(output_v(slidingVisX,slidingVisY,:), 3),kernel_av, 'same');
            input_va = zeros(visualFieldSize,visualFieldSize,visualLayers);
            input_va(slidingVisX,slidingVisY,:) = repmat(conv2( output_a,kernel_va, 'same'), [1, 1, visualLayers]) - c_va_gi_1_12 * sum(output_a(:));
            input_sa = conv2( output_a_fovSup, kernel_sa, 'same') - c_sa_gi_4 * sum(output_a(:));
            input_as = conv2(output_s, kernel_as, 'same');
            input_rg = c_rg_6_77 * output_g;
            input_rs = c_rs_6_78 * max(output_s(:));
            input_sr = c_sr_gi_5 * output_r;
            input_sx = c_sx_3_48 * output_x * (gauss2d_LAG1(spatialHalfSize,sigma_ss_exc_3_49,spatialHalfSize,spatialHalfSize,gazeShift(1),gazeShift(2)))';
            input_sg = c_sg_3_52 * output_g*output_s;
            input_ax = c_ax_2_20 * output_x *(gauss2d_LAG1(spatialHalfSize,sigma_aIn_2_22,spatialHalfSize,spatialHalfSize,gazeShift(1),gazeShift(2)))';
            input_ag = c_ag_2_28 * output_g * (gauss2d_LAG1(spatialHalfSize,sigma_aIn_2_22,spatialHalfSize,spatialHalfSize,gazeShift(1),gazeShift(2)))';
            input_ai = c_ai_2_21 * ior_input;% * output_i;
            input_xx = c_xx_5_71 * output_x;
            input_xfeature = c_xfeature_5_72*max(sigmoid(featureDetection_exc,beta_ft_x_73,inflec_ft_x_74));
            input_xr = c_xr_5_70 * output_r;
            input_gr = c_gr_4_64 * output_r;
            input_gx = c_gx_4_65 * output_x;
            input_rx = c_rx_6_76 * output_x;
            input_xg = c_xg_5_69 * output_g;
            input_fbButtonv = c_fv_7_82*fbButtonDetection;
            input_fbButtonfbButton = c_ff_7_83*output_fbButton;
            input_fbButtonExpfbButton = -output_fbButton;
            input_fbButtonExpfbButtonExp = -output_fbButtonExp;
            
            
            
            if fieldModel
                input_ff_exc = c_ff_7_83 * conv2( output_f,kernel_ff, 'same')' - 1;
                input_prefpref_exc = c_prefpref_exc_9_110 * conv2( output_f,kernel_prefpref, 'same')' - 1;
                input_cc_exc = c_cc_exc_8_99 * conv2( output_c,kernel_cc, 'same')';
            else
                input_ff_exc = c_ff_7_83 * output_f;
                input_prefpref_exc = c_prefpref_exc_9_110 * output_pref;
                input_cc_exc = c_cc_exc_8_99 * output_c;
            end
            
            
            % Feature Inihibition: "You will see this color or this color
            % in this specific location."
            if fieldModel
                input_ff_inh = 0;
                for i = 1:featureNum
                    input_ff_inh = input_ff_inh + -c_ff_inh_7_81*output_f(featureAntiCorrelations(i,2))*gauss([1:featureFieldSize],featureAntiCorrelations(i,4),sigma_ff);
                end
            else
                if size(featureAntiCorrelations) > 0
                    input_ff_inh(:,1) = -c_ff_inh_7_81*(output_f(featureAntiCorrelations(:,3)));
                    
                else
                    input_ff_inh = 0;
                end
                input_ff_inh = input_ff_inh - c_ff_gi_7_84*sum(output_f);
                input_ff_inh = input_ff_inh - ft_det_fbButton_gi*output_fbButton;
            end
            
            
            % Expectation Inhibition. If the model has seen a particular
            % feature value it should not want to see it again.
            if fieldModel
                input_preff = 0;
                maxFeatureActivations = max(output_f(featureAntiCorrelations(:,[2 4])'));
                for i = 1:featureNum
                    input_preff = input_preff + -c_preff_9_108*maxFeatureActivations(i)*gauss([1:featureFieldSize],output_f(featureAntiCorrelations(i,2)),sigma_ff);
                end
            else
                input_preff = kron(max(reshape(neurons_f,2,featureNum/2)),[1 1]);
                input_preff = -c_preff_9_108 * input_preff';
            end
            
            input_impPhasicChange = fbViewed * impatienceChangeStrength;
            
            input_impTimer = c_impTimer_11_123 * (trialTimeStep - (fbViewed*fbStartTime))^trial_impatience_exponent;
            
            neuron_imp = neuron_imp + (deltaT/tau)*(-neuron_imp + input_impTimer + input_impPhasicChange);
            
            input_gimp = c_gimp_4_67 * (trialTimeStep - fixationImp)^fixation_impatience_exponent;
            
            input_clickImp = c_clickimp_10_114 * output_imp;
            input_clickX = c_clickx_10_115 * output_x;
            
            input_clickc = c_clickc_10_113 * ((max(output_c(:))^2)/sum(output_c));
            %input_clickc = c_clickc_10_113 * max(output_c(:));
            
            input_gg = c_gg_4_66 * output_g;
            
            % This sums up all the inhibition from all the other categories.
            input_cc_inh(length(output_c),1) = 0;
            input_cc_exc(length(output_c),1) = 0;
            for input_counter = 1:length(output_c)
                %Could just sum and then subtract off what we don't want
                %in a future version.
                input_cc_inh(input_counter,1) = sum(-c_cc_inh_8_100*(output_c(~ismember(1:length(output_c),input_counter))));
                input_cc_exc(input_counter,1) = sum(c_cc_exc_8_99*(output_c(~ismember(1:length(output_c),input_counter))));
            end
            
            
            
            if fieldModel
                field_f = field_f +  (deltaT/tau) * (-field_f + h_f_80 + input_ff_exc' + input_ff_inh' + input_fv');
                field_pref = field_pref + (deltaT/tau) * (-field_pref + h_pref_106 + input_prefpref_exc' + input_prefc + input_preff');
                field_c = field_c(:) + (deltaT/tau) * (-field_c(:) + h_c_97 + input_cc_exc(:) + input_cc_inh + input_cf + feedback_c(:));
                
            else
                
                %% Feature Detection Neurons
                % Document reference: \ref{eq:ftdetectionsEq}
                % $$\tau \dot {f}(z,t) = -f(z,t) + {h_1}$$
                %
                neurons_f = neurons_f + (deltaT/tau) * (-neurons_f + h_f_80 + input_ff_exc + input_ff_inh + input_fv);
                
                %% Feature Expectation Neurons
                % Document reference: \ref{eq:ftexpNeuronsEq}
                % $$\tau \dot {f}(z,t) = -f(z,t) + {h_1}$$
                %
                
                neurons_pref = neurons_pref + (deltaT/tau) * (-neurons_pref + h_pref_106 + input_prefpref_exc + input_prefc + input_preff);
                
                %% Category Neurons
                % Document reference: \ref{eq:categoryNeuronsEq}
                % $$\tau \dot {c}(z,t) = -c(z,t) + {h_1}$$
                %
                neurons_c = neurons_c + (deltaT/tau) * (-neurons_c + h_c_97 + input_cc_exc + input_cc_inh + input_cf + feedback_c);
                
            end
            
            feedback_c = zeros(catFieldSize,1);
            
            
            %% Visual Field
            % Document reference: \ref{eq:visualField}
            % $$\tau \dot {v}(x,y,t) = -v(x,y,t) + {h_1} + r(x,t)  +  inp(x,y,t) + feature_{pref}(y,t) + \sum\nolimits_z {in{p_v}(x,y,z,t)} +  \int {\int {{w_v}(x - x',y - y',t)[{v^*}(x',y',t)]dx'dy'} } + {c_{12}}\sum\nolimits_{x,y} {{a^*}(x,y,t)} + {c_{13}} \sum\nolimits_z {{ft_{exp}}^*(x,y,z,t)} - {c_{14}}\sum\nolimits_z \int \int {{ft_{exp}}^*(x',y',z,t)})dx'dy' + {fb_{exp}}^*(t) + {\zeta _v}$$
            %
            % https://www.mathworks.com/matlabcentral/answers/161540-how-can-i-align-tex-equations-when-using-the-publish-functionality-of-the-matlab-editor
            %
            field_v = field_v + (deltaT/tau) * (-field_v + h_v_1 + input_v ...
                + input_vv + input_va + input_vpref);
            
            %% Spatial Attention Field
            % Document reference: \ref{eq:aField}
            % $$\tau \dot {a}(x,y,t) = - a(x,y,t) + {h_{a,18}} + \sum\nolimits_z {c_{19}{ft^{*}_{exp}}(z,t)v^{*}(x,y,z,t)} + c{}_{25}fb{_{\exp }^{*}}(t) + {c_{20}}{x^{*}}(t){\cal G}_A(fix({t_f}),{\sigma _{22}}){I_{sacc}} - {c_{21}}i(x,y,t) - {c_{28}}{g^{*}}(t){\cal G}_A(fix({t_f}),{\sigma _{22}}) - {c_{29}}{r^{*}}(t)a^{*}(x,y,t) + {c_{45}}{I_{sacc}}v_a(x,y,t) + s_a^{*}(x,y,t) + \zeta_{a}(x,y,t) + \int  \int  {w_a}(x - x',y - y',t){a^{*}}(x',y',t)dx'dy' - {c_{38}} {\int \int a^*(x',y',t) dx'dy'} $$
            %
            
            field_a = field_a + (deltaT/tau) * (-field_a + h_a_18 + input_ar  ...
                + input_aa + input_av + input_as + input_ax + input_ai + input_ag + attPref + input_aFBbutton);    %input_ai needs to get put back in.
            
            %% Saccade Motor Field
            % Document reference: \ref{eq:sField}
            % $$\tau \dot {v}(x,y,t) = - s(x,y,t) + {h_{s,47}} - {c_{48}}{x^*}(t){{\cal G}_S}(fix({t_f}),{\sigma _{49}}) + {c_{16}}{g^*}(t){s^*}(x,y,t) + a_s(x,y,t) + {\zeta _s}(x,y,t) + \int {\int {{w_s}(x - x',y - y',t){s^*}(x',y')dx'dy'} } - {c_{62}}\sum\nolimits {{s^*}(x,y,t)}$$
            %
            field_s = field_s + (deltaT/tau) * (-field_s + h_s_47 ...
                + input_ss + input_sa + input_sr + input_sx + input_sg);
            
            
            %% Gaze Change Neuron
            % Document reference: \ref{eq:gazeChangeEq}
            %
            neuron_g = neuron_g + (deltaT/tau) * (-neuron_g + h_g_63 + input_gr + input_gg + input_gimp + input_gx);
            
            %% Saccade Fixation Neuron
            % Document reference: \ref{eq:fixNeuronEq}
            %
            neuron_x = neuron_x + (deltaT/tau) * (-neuron_x + h_x_68 + input_xr + input_xg + input_xx + input_xfeature); %input_xChanges
            
            
            %% Saccade Reset/Initiation Neuron
            % Document reference: \ref{eq:saccInitiationEq}
            %
            neuron_r = neuron_r + (deltaT/tau) * (-neuron_r + h_r_75 + input_rr + input_rs + input_rg + input_rx);
            
            
            %% Click Decision Neuron
            % Document reference: \ref{eq:visualField}
            %
            neuron_click = neuron_click + (deltaT/tau) * (output_click - neuron_click + input_clickc +input_clickImp + input_clickX);
            
            
            
            
            neuron_fbButton = neuron_fbButton + (deltaT/tau) * (-neuron_fbButton + input_fbButtonv + input_fbButtonfbButton);
            neuron_fbButtonExp = neuron_fbButtonExp + (deltaT/tau) * (-neuron_fbButtonExp + input_fbButtonExpfbButton + input_fbButtonExpfbButtonExp);
            
            if logging
                GazeLvl = [GazeLvl;subjectNumber TETTimeTimeSteps TETTime gazeShift 0 0 trialNum trialPhase];
            end
            
            TETTime = TETTime +  deltaT ; % Eyetracker timestamp units. Timesteps from start of experiment.
            TETTimeTimeSteps = TETTimeTimeSteps + 1;
            

            if noise
                if fieldModel
                    field_f = field_f + 0* ((sqrt(deltaT)/sqrt(tau))*q_f_95) *randn(featureFieldSize,1);
                    field_pref = field_pref + ((sqrt(deltaT)/sqrt(tau))*q_pref_111) *randn(featureFieldSize,1);
                    field_c = field_c + 0* ((sqrt(deltaT)/sqrt(tau))*q_c_103) .*randn(catFieldSize,1);
                else
                    neurons_f = neurons_f + ((sqrt(deltaT)/sqrt(tau))*q_f_95) *randn(featureNum,1)*q_spatial_f_96;
                    neurons_pref = neurons_pref + ((sqrt(deltaT)/sqrt(tau))*q_pref_111) *randn(featureNum,1)*q_spatial_pref_112;
                    neurons_c = neurons_c + ((sqrt(deltaT)/sqrt(tau))*q_c_103) .*randn(catNum,1)*q_spatial_c_104;
                end
                for i=1:visualLayers
                    % sqrt(gridRefinement)^2 is needed to compensate for
                    % grid spacing effects and power of 2 is because we're
                    % in 2 dimensions.
                    field_v(:,:,i) = field_v(:,:,i) + conv2( ((sqrt(deltaT)/sqrt(tau))*q_v_5) * sqrt(gridRefinement)^2 * randn(visualFieldSize, visualFieldSize), kernel_q_spatial, 'same');
                end
                field_a = field_a + conv2( ((sqrt(deltaT)/sqrt(tau))*q_a_39) * sqrt(gridRefinement)^2 * randn(spatialFieldSize, spatialFieldSize),kernel_q_spatial, 'same')*q_spatial_a_40;
                field_s = field_s + conv2( ((sqrt(deltaT)/sqrt(tau))*q_s_57) * sqrt(gridRefinement)^2 * randn(spatialFieldSize, spatialFieldSize), kernel_q_spatial, 'same')*q_spatial_s_58;
                neuron_r = neuron_r + ((sqrt(deltaT)/sqrt(tau))*q_r_93) * randn*q_spatial_r_94;
                neuron_x = neuron_x + ((sqrt(deltaT)/sqrt(tau))*q_x_91) * randn*q_spatial_x_92;
                neuron_g = neuron_g + ((sqrt(deltaT)/sqrt(tau))*q_g_89) * randn*q_spatial_g_90;
                neuron_imp = neuron_imp + ((sqrt(deltaT)/sqrt(tau))*q_imp_124) * randn*q_spatial_imp_125;
                neuron_click = neuron_click + ((sqrt(deltaT)/sqrt(tau))*q_click_118) *randn*q_spatial_click_119;
            end
            
            if trialNum == 1 || trialNum == 130 || cpDetector
                neuronHistory = [neuronHistory; trialNum trialPhase output_g output_x output_r output_f' output_fbButton trialFix length(trialFixations) max(input_fv) output_c'];
            end
            
            if trialNum==181
                save('neuronHistory.mat','neuronHistory')
            end
            
            if visualize == 1
                
                if fieldModel
                    set(hPlot_f, 'XData', output_f);
                    set(hPlot_pref, 'XData', output_pref);
                    set(hPlot_c, 'XData', output_c);
                else
                    set(Ft_Det_Axis, 'YData', output_f);
                    set(Ft_Exp_Axis, 'YData', output_pref);
                    set(Category_Neuron_Axis, 'YData', output_c);
                end
                % three 2d fields are plotted with transpose so that x is horiz
                % and y is vert.
                
                set(Visual_Field_Axis, 'CData', flipud(sum(output_v,3)'));
                set(Attention_Field_Axis, 'CData', output_a');
                set(Saccade_Field_Axis, 'CData', output_s');
                if ~environLoaded
                    set(Spatial_Gain_Axis, 'CDdata', attPref');
                end
                if maxFields
                    axes(Max_Att_Axis)
                    cla
                    axes(Max_Sacc_Axis)
                    cla
                    axes(Raw_Max_Sacc_Axis)
                    cla
                    axes(Raw_Max_Att_Axis)
                    cla
                    
                    line([1:51],max((field_s')),'Parent',Raw_Max_Sacc_Axis,'Color','k')
                    line([1:51],max((field_a')),'Parent',Raw_Max_Att_Axis,'Color','k')
                    
                    line(Max_Sacc_Axis,[1:51],max(fliplr(output_s')),'Color','b')
                    line(Max_Att_Axis,[1:51],max(fliplr(output_a')),'Color','b')
                    %
                    
                end
                
                set(Weight_Axis, 'CData', output_wT');
                set(F1_circle, 'XData', [stimPos_x(1)], 'YData', [stimPos_y(1)], 'MarkerFaceColor', featureColours(:,1)', 'MarkerEdgeColor', featureColours(:,1)', 'MarkerSize', 10);hold on
                if size(featureColours,2) >1
                    set(F2_circle, 'XData', [stimPos_x(2)], 'YData', [stimPos_y(2)], 'MarkerFaceColor', featureColours(:,2)', 'MarkerEdgeColor', featureColours(:,2)', 'MarkerSize', 10);hold on
                end
                if size(featureColours,2) >2
                    set(F3_circle, 'XData', [stimPos_x(3)], 'YData', [stimPos_y(3)], 'MarkerFaceColor', featureColours(:,3)', 'MarkerEdgeColor', featureColours(:,3)', 'MarkerSize', 10);hold on
                end
                %In case we're visualizing 4 features.
                if locationNum == 4
                    set(Screen_Display_Axis(6),'XData', [stimPos_x(4)], 'YData', [stimPos_y(4)], 'MarkerFaceColor', featureColours(:,4)', 'MarkerEdgeColor', featureColours(:,4)', 'MarkerSize', 10);hold on
                    %set(Screen_Display_Axis(7), 'YData', [stimPos_y(4)], 'MarkerFaceColor', featureColours(:,4)', 'MarkerEdgeColor', featureColours(:,4)', 'MarkerSize', 10);hold on
                else
                    set(Screen_Display_Axis(6), 'XData', [spatialFieldFovea(1)-spatialHalfSize], 'YData', [spatialFieldFovea(2)-spatialHalfSize] , 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 20);
                end
                
                set(Gain_Applied_Axis, 'YData', c_prefgain_9_105*featureValue'); %input_apref
                
                
                set(Trial_Imp_Axis, 'YData', output_imp);
                set(Motor_Click_Axis, 'YData', output_click);
                set(Saccade_Neuron_Axis, 'YData', output_r);
                set(Fixation_Neuron_Axis, 'YData', output_x);
                set(Gaze_Change_Axis, 'YData', output_g);
                set(FB_Det_Axis, 'YData', output_fbButton);
                set(FB_Exp_Axis, 'YData', output_fbButtonExp);
                
                set(Trial_Num_Axis,'String', ['Trial = ' num2str(trialNum)])
                set(Exp_Time_Axis,'String', ['Experiment time (ms) = ' num2str(round(TETTime))])
                
                pause(0.02);
            end
            
            %If the end of the trial is detected...
            if trialTimeStep == tMax
                
                trialFixationsTemp = trialFixations;
                trialFixationsTemp(trialFixations==1)=LocationRelevance.one;
                trialFixationsTemp(trialFixations==2)=LocationRelevance.two;
                
                if featureNum ==6
                    trialFixationsTemp(trialFixations==3)=LocationRelevance.three;
                elseif featureNum ==8
                    trialFixationsTemp(trialFixations==4)=LocationRelevance.four;
                end
                trialFixations = trialFixationsTemp;
                
                if ColourOnset ==0
                    ColourOnset = TETTime;
                end
                
                if firstFBTrialFlag
                    
                    
                    trialLogger{2,14} = [trialLogger{2,14}; CorrectResponse{end}];
                    trialLogger{2,15} = [trialLogger{2,15}; Response];
                    trialLogger{2,18} = [trialLogger{2,18}; TrialAccuracy];
                    %Log the category and feature activations
                    trialLogger(2,1:4) = [{[trialLogger{2,1} neurons_c]} {[trialLogger{2,2} output_c]} {[trialLogger{2,3} neurons_f]} {[trialLogger{2,4} output_f]}];
                    trialLogger(2,5:6) = [{[trialLogger{2,15}; {wT}]} {[trialLogger{2,16}; {output_wT}]}];
                end
                
                FeedbackRT = TETTime - ColourOnset;
                FeedbackOnset = ColourOnset;
                
                
                
                %This is to signify that transfer should begin based on 18
                %trials in a row being correct.
                
                if cpDetector == 0 && ~feedbackOff
                    if (trialNum > 18) && (sum(trialLogger{2,18}(trialNum-18:trialNum-1,1)) == 18)
                        cpDist = NaN;
                        cp = trialNum;
                        cpDetector = 1;
                    end
                end
                
                
                if logging
                    if locationNum ==3
                        TrialLvl = [subjectNumber trialNum binaryFeatureValues(LocationRelevance.one) binaryFeatureValues(LocationRelevance.two)	binaryFeatureValues(LocationRelevance.three) CorrectResponse{end} Response	TrialAccuracy	StimulusRT	FeedbackRT	1	FixationOnset	StimulusOnset	ColourOnset	fixCrossInd	fixCrossInd	FeedbackOnset];
                    elseif locationNum ==4   
                        TrialLvl = [subjectNumber trialNum binaryFeatureValues(LocationRelevance.one) binaryFeatureValues(LocationRelevance.two)	binaryFeatureValues(LocationRelevance.three) binaryFeatureValues(LocationRelevance.four) CorrectResponse{end} Response	TrialAccuracy	StimulusRT	FeedbackRT	1	FixationOnset	StimulusOnset	ColourOnset	fixCrossInd	fixCrossInd	FeedbackOnset];
                                                
                    end
                    dlmwrite(['./dft' expName 'TrialLvl-' num2str(subjectNumber) '.txt'], TrialLvl,'delimiter', '\t', 'precision','%.0f', '-append');
                end
            end
            
            
            %If the end of the experiment is detected...
            if trialTimeStep == tMax && trialNum == maxTrials
                
                trialLogger{2,16} = [trialLogger{2,16};trialTimeStep];
                trialLogger{2,10} = [trialLogger{2,10}; {trialFixations'}];
                if ~feedbackOff
                    trialLogger{2,11} = [trialLogger{2,11}; {deltaT*(saccStartTimes(2:length(trialFixations)+1) - saccEndTimes(1:length(trialFixations)))}];
                    trialLogger{2,12} = [trialLogger{2,12}; {trialFixationsFB'}];
                    trialLogger{2,13} = [trialLogger{2,13}; {deltaT*(saccStartTimes(length(trialFixations)+2:end) - saccEndTimes(length(trialFixations)+1:end-1))}];
                    trialLogger{2,17} = [trialLogger{2,17}; {deltaT*(saccEndTimes - saccStartTimes)}];
                else
                    trialLogger{2,11} = [trialLogger{2,11}; {deltaT*(saccStartTimes(2:length(trialFixations)+1) - saccEndTimes(1:length(trialFixations)))}];
                    trialLogger{2,12} = [trialLogger{2,12}; {nan}];
                    trialLogger{2,13} = [trialLogger{2,13}; {nan}];
                    trialLogger{2,17} = [trialLogger{2,17}; {deltaT*(saccEndTimes - saccStartTimes)}];
                end
                trialLogger{2,19} = [trialLogger{2,19}; {accuracyLevels}];
                trialLogger{2,21} = [trialLogger{2,21}; {fixationChangesRecord}];
                trialLog(trialLogger);
                
                if structure == 1
                    
                    if (category(trialNum) == 1) || (category(trialNum) ==2) && structure ==1
                        display('Did the model look at the irrelevant feature?')
                        ismember(LocationRelevance.three,trialFixations)
                        trialLogger{2,20} = [trialLogger{2,20}; ismember(LocationRelevance.three,trialFixations)];
                    else
                        display('Did the model look at the irrelevant feature?')
                        ismember(LocationRelevance.two,trialFixations)
                        trialLogger{2,20} = [trialLogger{2,20}; ismember(LocationRelevance.two,trialFixations)];
                    end
                end
                
                category = [trainingCategories;category];
                
                if visualize
                    stopPlay;
                end
                
                if logging
                    if locationNum ==3
                        ExpLvl = [ExpLvl; subjectNumber feedbackType LocationFeature.one LocationFeature.two LocationFeature.three LocationRelevance.one LocationRelevance.two LocationRelevance.three];
                    elseif locationNum ==4
                        ExpLvl = [ExpLvl; subjectNumber feedbackType LocationFeature.one LocationFeature.two LocationFeature.three LocationFeature.four LocationRelevance.one LocationRelevance.two LocationRelevance.three LocationRelevance.four];
                    end
                    
                    dlmwrite(['./dft' expName 'ExpLvl.txt'], ExpLvl,'delimiter', '\t', 'precision','%.12g', '-append');
                    x = genvarname('trialLogger');
                    eval(['save', ' trialLog' num2str(subjectNumber) '.mat ', x])
                    dlmwrite(['./dft' expName 'GazeLvl-' num2str(subjectNumber) '-1.gazedata'], GazeLvl,'delimiter', '\t', 'precision','%.0f', '-append');
                end
                initExperiment = false;
                return;
                
            elseif (trialTimeStep == tMax) && trialNum <= maxTrials
                trialLogger{2,16} = [trialLogger{2,16};trialTimeStep];
                trialLogger{2,10} = [trialLogger{2,10}; {trialFixations'}];
                if ~feedbackOff
                    trialLogger{2,11} = [trialLogger{2,11}; {deltaT*(saccStartTimes(2:length(trialFixations)+1) - saccEndTimes(1:length(trialFixations)))}];
                    trialLogger{2,12} = [trialLogger{2,12}; {trialFixationsFB'}];
                    trialLogger{2,13} = [trialLogger{2,13}; {deltaT*(saccStartTimes(length(trialFixations)+2:end) - saccEndTimes(length(trialFixations)+1:end-1))}];
                    trialLogger{2,17} = [trialLogger{2,17}; {deltaT*(saccEndTimes - saccStartTimes)}];
                else
                    trialLogger{2,11} = [trialLogger{2,11}; {deltaT*(saccStartTimes(2:length(trialFixations)+1) - saccEndTimes(1:length(trialFixations)))}];
                    trialLogger{2,12} = [trialLogger{2,12}; {nan}];
                    trialLogger{2,13} = [trialLogger{2,13}; {nan}];
                    trialLogger{2,17} = [trialLogger{2,17}; {deltaT*(saccEndTimes - saccStartTimes)}];
                end
                trialLogger{2,21} = [trialLogger{2,21}; {fixationChangesRecord}];
                trialLog(trialLogger);
                
                if structure == 1
                    if (category(trialNum) == 1) || (category(trialNum) ==2) && structure ==1
                        display('Did the model look at the irrelevant feature?')
                        ismember(LocationRelevance.three,trialFixations)
                        trialLogger{2,20} = [trialLogger{2,20}; ismember(LocationRelevance.three,trialFixations)];
                    else
                        display('Did the model look at the irrelevant feature?')
                        ismember(LocationRelevance.two,trialFixations)
                        trialLogger{2,20} = [trialLogger{2,20}; ismember(LocationRelevance.two,trialFixations)];
                    end
                end
                
                initTrial = 1;
                trialTimeStep = trialTimeStep + 1;
                toc
            end
            trialTimeStep = trialTimeStep + 1;
            
            if ~environLoaded
                load Trial359state.mat
                environLoaded = 1;
                
                hAllAxes = findobj(gcf,'type','axes')
                
                Visual_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Visual Field'),get(hAllAxes,'Title')))).Children;
                Attention_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Spatial Attention Field'),get(hAllAxes,'Title')))).Children;
                Saccade_Field_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Saccade Motor Field'),get(hAllAxes,'Title')))).Children;
                
                hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Visual Field'),get(hAllAxes,'Title')))).CLim=[0 1];
                hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Spatial Attention Field'),get(hAllAxes,'Title')))).CLim=[0 1];
                hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Saccade Motor Field'),get(hAllAxes,'Title')))).CLim=[0 1];
                Gaze_Change_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Gaze Change'),get(hAllAxes,'Title')))).Children;
                Fixation_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Fixation Excitation'),get(hAllAxes,'Title')))).Children;
                Saccade_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Saccade Reset'),get(hAllAxes,'Title')))).Children;
                FB_Exp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'FB  Expectation'),get(hAllAxes,'Title')))).Children;
                FB_Det_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'FB detection'),get(hAllAxes,'Title')))).Children;
                Category_Neuron_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Category Neurons'),get(hAllAxes,'Title')))).Children;
                Ft_Exp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Feature Expectations'),get(hAllAxes,'Title')))).Children;
                Ft_Det_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Feature Detection'),get(hAllAxes,'Title')))).Children;
                Motor_Click_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Motor Click Neuron'),get(hAllAxes,'Title')))).Children;
                Trial_Imp_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Trial Impatience'),get(hAllAxes,'Title')))).Children;
                Trial_Num_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children(1);
                Exp_Time_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children(2);
                F1_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children(3);
                F2_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children(4);
                F3_circle = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children(5);
                Max_Att_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Max Attention Activity'),get(hAllAxes,'Title'))));
                Weight_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Weight Matrix'),get(hAllAxes,'Title')))).Children;
                Max_Sacc_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Max Saccade Activity'),get(hAllAxes,'Title'))));
                Screen_Display_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Screen w/fovea'),get(hAllAxes,'Title')))).Children;
                
                Raw_Max_Att_Position = Max_Att_Axis.Position; % position of first axes
                Raw_Max_Att_Axis = axes('Position',Raw_Max_Att_Position,'XAxisLocation','top',...
                    'YAxisLocation','left',...
                    'Color','none','ylim',[-5.5 0],'YTick',[],'xlim',[1 51],'XTick',[]);
                
                Raw_Max_Sacc_Position = Max_Sacc_Axis.Position; % position of first axes
                Raw_Max_Sacc_Axis = axes('Position',Raw_Max_Sacc_Position,...
                    'XAxisLocation','top',...
                    'YAxisLocation','left',...
                    'Color','none','ylim',[-6 13],'YTick',[],'xlim',[1 51],'XTick',[]);
                
                Gain_Applied_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,...
                    'Gain'),get(hAllAxes,'Title')))).Children;
                
                Spatial_Gain_Axis = hAllAxes(find(cellfun(@(x) strcmp(x.String,'Space Coded'),get(hAllAxes,'Title')))).Children;
                
                
                
                
            end
            
        end
        
        if trialNum ==1 && profiling
            profile viewer
            p = profile('info');
            profsave(p,'profile_results')
        end
        
        if logging
            dlmwrite(['./dft' expName 'GazeLvl-' num2str(subjectNumber) '-1.gazedata'], GazeLvl,'delimiter', '\t', 'precision','%.0f', '-append');
            GazeLvl=[];
        end
        trialNum = trialNum + 1;
        
        %Setting options for the the switch to the transfer phase.
        if (transfer.switch && (trialNum ==maxTrials) && ~(feedbackOff)) || (transfer.switch && cpDetector && ~feedbackOff)
            trialNum = 1001;
            feedbackOff = 1;
            trainingCategories = category;
            category = transfer.categories;
            stimulus = transfer.stims;
            transfer.switch = 2;
            maxTrials = 1032;
            maxTrials = 1004;
        end
        
        pause(0.02);
        
    end
    pause(0.02);
    
end


    function playButtonCallback(hObject, event)
        if state == PAUSE
            %dbquit;
            %eval('return');
            startPlay
        elseif state == PLAY
            pausePlay;
        elseif state == STOP
            initExperiment = true;
            initTrial = true;
            startPlay;
        end
    end

    function stopButtonCallback(hObject, event)
        stopPlay;
    end


    function startPlay
        state = PLAY;
        %set(playButton, 'CData', pauseIcon);
    end

    function pausePlay
        state = PAUSE;
        %set(playButton, 'CData', playIcon);
        keyboard
    end

    function stopPlay
        state = STOP;
        %set(playButton, 'CData', playIcon);
    end


    function figureCloseRequest(hObject, event)
        state = QUIT;
        closereq;
    end


end