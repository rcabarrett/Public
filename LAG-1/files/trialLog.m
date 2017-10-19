%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 0r0v
%   - this should probably be replaced with a standard modeling display wrapper 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trialLog(metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Trial Statistics Display %%
%
% Author: Jordan
% Date: Sep, 2013
% Brief description:
%   Simply displays experiment information on a trial by trial basis.
%
% 
% trialLogger(1:19) = [{'Category neurons at decision time'},...
%     {'Sigmoided Category neurons at decision time'}, {'Feature neurons at decision time'},...
%     {'Sigmoided Feature neurons at decision time'}, {'Weights at decision time'},...
%     {'Sigmoided weights at decision time'}, {'Phase 2 Reaction time'},...
%     {'Phase 4 Reaction time'}, {'Raw weight changes'},...
%     {'Phase 2 Fixations'}, {'Phase 2 Fixation durations'},...
%     {'Phase 4 Fixations'}, {'Phase 4 Fixation durations'}, {'Correct Category'},...
%     {'Response'}, {'Total Trial Time'}, {'Saccade Times'},{'Accuracy'}, {'Accuracy Levels'}];
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Review: X
% Verify: Y





display(['The right answer was ' num2str(metrics{2,14}(end)) ]) % display the right answer.
display(['The response was ' num2str(metrics{2,15}(end)) ]) % display the model's response.

display('The category activations at time of decision selection were:')
printmat([metrics{2,1}(:,end) metrics{2,2}(:,end)])

display('The weights at time of decision selection were:')
printmat(metrics{2,5}{end})

display('The phase 2 fixation locations and durations were:')
printmat([metrics{2,10}{end}';metrics{2,11}{end}'])

display('The phase 4 fixation locations and durations were:')
printmat([metrics{2,12}{end}';metrics{2,13}{end}'])

display('The fixation changes were:')
printmat([metrics{2,21}{end}])

display('The accuracy was:')
printmat([metrics{2,18}(end)])
end
