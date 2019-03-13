function dC = Css(t,C)
% Steady state serum concentration 2-compartment model
% find Css for 3 ODEs representing the compartments, the first 
% compartment, the GI tract, is trivial.

dC = zeros(3,1);

% Rate constants
k_a = 14;
k10 = 2.9;
k12 = 9.5;
k21 = 2.1;
V1 = 0.9;

dC(1) = -k_a*C(1); % dA_gi/dt
dC(2) = (k_a*C(1))/V1 - (k10+k12)*C(2) + k21*C(3); % dC1/dt
dC(3) = k12*C(2) - k21*C(3); % dC2/dt