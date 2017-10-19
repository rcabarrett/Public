%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 9a3b2c0r0v
%   - Input/Output data structure decimation needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wT = associator1(wT, tau, deltaT, pre_synaptic_neurons, post_synaptic_neurons, learning_rate, TrialAccuracy, feedbackType,activation_threshold_presynaptic,activation_threshold_postsynaptic,antiHebbianFlag, decay_rate, antiHebbRate, CategoryToAssociate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% associator.m %%
%
% Author: Calen Walshe, Jordan B.
% Date: Mar, 2013
%
%   Brief description:
%       Updates a weight matrix one time, based on the learning parameters
%       supplied.
%
%
% Inputs:
% 
% wT = the weight matrix to be updated.
% tau = the time scale of the 
% pre_synaptic_neurons = the activation of the pre-synaptic neurons
% post_synaptic_neurons = the activation of the post-synaptic neurons
% learning_rate = the desired level of association per time step.
% TrialAccuracy = needed if error is going to entail anything.
% feedbackType = 
% 
% % From LAG-1:
% 
% %   feedbackType = 1 = corrective feedback
% %   feedbackType = 2 = learning only on correct trials
% %   feedbackType = 3 = anti-Hebbian learning on error trials
% %   feedbackType = 4 = anti-Hebbian and Hebbian learning on error trials
% 
% activation_threshold_presynaptic = this should be ammended to the sigmoided
% output level in a future version, instead of this arbitrary threshold for
% participation in the associative procedure.
% activation_threshold_postsynaptic = ditto for this.
% antiHebbianFlag = turn on anti-Hebbian association
% decay_rate = weakens connections at a constant rate
% deltaT = the characteristic time of the system
% antiHebbRate = the desired level of de-association per time step
% CategoryToAssociate = the specific category undergoing association
% 
% Outputs:
% 
% wT = the updated weight matrix.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ErrorBoost = 10;

postMat=repmat(post_synaptic_neurons,1,length(pre_synaptic_neurons));
preMat=repmat(pre_synaptic_neurons',length(post_synaptic_neurons),1);

        % positive association
  postI=post_synaptic_neurons>activation_threshold_postsynaptic;
  preI=pre_synaptic_neurons>activation_threshold_presynaptic;
%  wT(postI, preI) = wT(postI, preI) ... 
%                + (learning_rate*(deltaT/tau)) * ...
%                ( -wT(postI, preI).* preMat(postI,preI) + postMat(postI,preI).*preMat(postI,preI));
            
  wT(postI, preI) = wT(postI, preI) ... 
                + (learning_rate*(deltaT/tau)) * ...
                (postMat(postI,preI).*preMat(postI,preI));
              
        
        % weight decay
        %wT = wT + learning_rate/decay_rate*deltaT/tau *wT; %growth
        wT = wT - learning_rate/decay_rate*deltaT/tau *wT; %decay
        
        
        % additional anti-Hebbian learning for the
        % corrective feedbackType on error trials.
        if TrialAccuracy==0 && antiHebbianFlag ==1
           
            psn1=ErrorBoost*ones(length(post_synaptic_neurons),1);
            psn1(CategoryToAssociate,1)=-ErrorBoost;
            %postI1=post_synaptic_neurons<activation_threshold_postsynaptic;
            postI1=psn1<activation_threshold_postsynaptic; %Jb change Oct 3,2016
                      
                wT(postI1, preI) = wT(postI1, preI)...
                    + (-(learning_rate/antiHebbRate)*(deltaT/tau)) * (-wT(postI1, preI) ...
                + abs(kron(psn1(postI1),pre_synaptic_neurons(preI)')));
            
        end
            
        
    end
    
   
