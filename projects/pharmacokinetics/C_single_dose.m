function dC = C_single_dose(t,C,flag,SA,k_a,k10,k12,k21,V1)
% Single dose dynamics of serum concentration 2-compartment model
% find Css for 3 ODEs representing the compartments, the first 
% compartment, the GI tract, is trivial.

dC = zeros(3,1);
if SA
  % dont assign  parameters if we are in Sensitivity Analysis mode, they 
  % will be assigned during function call from SensitivityAnalysis.m
else
    % Otherwise, assign parameters according to reference [7] Table 1
    k_a = 14;  % absorption
    k10 = 2.9; % elimination
    k12 = 9.5; % distribution
    k21 = 1.4; % redistribution
    V1 = 0.9; % volume
end 

dC(1) = -k_a*C(1); % dA_gi/dt   iv dosage
dC(2) = (k_a*C(1))/V1 - (k10+k12)*C(2) + k21*C(3); % dC1/dt  oral dosage
dC(3) = k12*C(2) - k21*C(3); % dC2/dt