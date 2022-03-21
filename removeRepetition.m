% *******************************************************
% Function removeRepetition()
% removeRepetition() deletes the repeated values in the 
% time and states vectors created in the segmentation 
% process (while enforcing continuity conditions)
% Ouput: [T,X]
% T: time vector without repeated values
% X: states struct without repeated values
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function [T,X] = removeRepetition(X,te)
T               = te(:);
[~, ind]        = unique(T);
duplicate_ind   = flip(setdiff(1:size(T), ind));
duplicate_value = T(duplicate_ind);
for i = 1 : length(duplicate_ind)
    T(duplicate_ind(i))    = [];
    X.X1(duplicate_ind(i)) = [];
    X.X2(duplicate_ind(i)) = [];
    X.L1(duplicate_ind(i)) = [];
    X.L2(duplicate_ind(i)) = [];
end
assignin('base','repetitions_in_sol',duplicate_value)
assignin('base','duplicate_ind',duplicate_value')