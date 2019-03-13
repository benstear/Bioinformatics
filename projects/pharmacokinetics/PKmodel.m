%
%
% Compartmental Kinetics

% Kinetic rate constants
% Dosage of 1.1 (mg/kg) 
k_a = 14; % absorption
k10 = 3.3; % same as elimination, k_e
k12 = 9.5;
k21 = 2.1;
V1 = 0.9;

% Dose (mg)
D = 20;
% Salt factor
S = 0.92;
% Bio-availability, for IV administration, F=1.0,
F = 0.5; 
% Volume of distribution (L/kg)
Vd = 8.7;

% initial conditions
Agi_t0 = D*S*F; % Amount of drug in GI track at time = 0 (effective dose)
C1_t0= 0;
C2_t0 = 0;
C0 = [Agi_t0 C1_t0 C2_t0];
tspan = 1:0.01:10; 

% Single dose model: solve system and plot concentration vs time for each
% compartment SA,...
SA = 0; % we are not in SA mode
                            
[t, C] = ode45('C_single_dose',tspan,C0,SA,k_a,k10,k12,k21,V1); 
%{                                          
plot(t,C(:,1),'LineWidth',3)  % ******** GI tract ******* 
hold on
plot(t,C(:,2),'LineWidth',3)  % ********Compartment 1 *****
                              % (central compartment, 
                              % represents serum and tissue concentration,
                              % where absorption is basically instantaneous)

plot(t,C(:,3),'LineWidth',3)  % ****** Compartment 2 (peripheral compartment,
                               % represents tissues where absorption of the
                               % drug is slower)

lgnd = legend('GI tract','Serum','Peripheral Tissue'); lgnd.FontSize = 16;
%title('Pharmacokinetic profile of Paroxetine (Single dose dynamics)')
title('Two Compartment Simulation of Paroxetine (Single dose)')
ylim([0 10]), xlim([1 10]) %;axis('tight')
hold off;
xlabel('Time (hours)'), ylabel('Concentration (ng/mL)')
%}


% Calculate AUC 
AUC_iv = trapz(t,C(:,1)); % GI tract
AUC_oral = trapz(t,C(:,2)); % serum

% Calculate half-life 
Cmax_iv = max(C(:,1));
half_life_val = Cmax_iv/2; % 3.022
t_half_iv = t(12); %------- 1.11

Cmax_oral = max(C(:,2));
half_life_val = Cmax_oral/2; % 3.08
t_half_oral = t(25); % -----1.24

% Calculate Tmax
%Tmax_oral = %t(Cmax_oral);

% Create table of parameters comparing predicted values and experimental
% data (For mice)
params = {'AUC';'C_max';'t_1/2'};
Predicted_IV = [AUC_iv;Cmax_iv;t_half_iv];
Experimental_IV = [];
Predicted_oral = [AUC_oral;Cmax_oral;t_half_oral];
Experimental_oral = [];
T = table(params,Predicted_IV,Predicted_oral);


%% Validation with standard pharmacokinetic equations
C_single_val = zeros(1,length(tspan));
for i = 1:length(tspan)
    % Compute GI tract serum concentration of one dose.
    C_single_val(i) = ((F*D*k_a)/Vd*(k_a-k10))*(exp(-k10*tspan(i))-exp(-k_a*tspan(i)));
end
plot(t,C(:,1),'LineWidth',3), hold on
plot(tspan,C_single_val,'LineWidth',3)
lgnd = legend('Predicted','Experimental'); lgnd.FontSize = 14;
title('Prediction vs Experimental serum concentration (IV)')
ylim([0 6]), xlim([1 6]);
xlabel('Time (h)'), ylabel('Concentration (ng/M)')

mse = immse(C(:,1),C_single_val');


%t_max_single_dose = log(k_a/k10)/(k_a - k10);
%plot(1,t_max_single_dose,'*k')

%% Steady state modeling

% The mean elimination half-life is approximately 21 hours (human)
% after oral dosing of Paroxetine tablets, 30 mg daily for 30 days.

%{
dosing_interval = 24; % once daily
tau = D/dosing_interval;
C_multiple_dose = [];

for i = 1:length(tspan)
    C_multiple_dose(i) = Agi_t0*exp(-k10*tspan(i))/(1-exp(-k10*tau));
end
plot(tspan,C_multiple_dose,'LineWidth',3)

%}
t_half = log(2)/k10;




%{
https://reference.medscape.com/drug/paxil-brisdelle-paroxetine-342959#0
20-30mg is a typical daily dose
Absorption:
Peak plasma time: 5.2-8.1 hr (immediate-release); 6-10 hr (controlled-release)
Peak plasma concentration: 61.7 ng/mL
Distribution:
Protein bound: 93-95%
Vd: 8.7 L/kg

WIKI:
 It has an absolute bioavailability of about 50%, with evidence of a 
 saturable first-pass effect.[70] When taken orally, it achieves maximum 
 concentration in about 6?10 hours[63] and reaches steady-state in 7?14 
 days
%}