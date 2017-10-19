%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Revision Code: 1r0v
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Softer max is an implementation of the softmax decision rule that employs
% a temperature parameter for use in reinforcement learning.

% Takes temperature, possibleActionVect as input; outputs the probability
% of selecting an action.

%Author: Caitlyn McColeman
%Date: February 21 2011

%Reviewed: Jordan B.
%Verified: Paul (changed the location of the ) and advising to produce
%normalized probability vector.

function P = softerMax(possibleActionVect, temp)

for action=1:length(possibleActionVect)
    P(action) = (exp(possibleActionVect(action)/temp))./sum(exp(possibleActionVect/temp));
end