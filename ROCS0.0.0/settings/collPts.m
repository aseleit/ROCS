% *******************************************************
% Function collPts()
% Sets the collocation parameters
% N:        number of nodes
% M:        number of segments
% Neq:      number of equations to be solved
% tf:       problem final time [sec]
% Output:   struct gridparam with fields
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function gridparam = collPts()
gridparam.N         = 40;
gridparam.M         = 3; 
gridparam.Neq       = 4;
gridparam.tf        = 2;
gridparam.RBFtype   = 'MCQ';
end