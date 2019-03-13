
%% Sensitivity Analysis
% Order by importance the strength and relevance of the inputs 
% in determining the variation in the output.

SA=1; % Indicate we want to run the model in sensitivity analysis mode.

% The range of values that we will test the parameters on
% will be: k_n = [n*0.1 n/4 n/2 n n*2 n*4 n*10]
% Kinetic rate constants
% Dosage of 1.1 (mg/kg) 
k_a = 14; % absorption
k10 = 3.3; % same as elimination, k_e
k12 = 9.5;
k21 = 2.1;
V1 = 0.9;
tspan = 1:0.01:10; 

D = 20;
S = 0.92;
F = 0.5; 
Agi_t0 = D*S*F;
C1_t0= 0;
C2_t0 = 0;
C0 = [Agi_t0 C1_t0 C2_t0];

SA_range = [0.1 1/4 1/2 1 2 4 10];
k_a_range = k_a*SA_range; 
k10_range = k10*SA_range;
k12_range = k12*SA_range;
k21_range = k21*SA_range;
V1_range = V1*SA_range;
param_ranges = [k_a_range; k10_range; k12_range; k21_range; V1_range];

title('Sensitivity Analysis for the 4 kinetic rates')
% Plot C(:,1), the GI tract concentration/IV administration
for i = 1:5
    for j = 1:7
    
    if i == 1 % vary k_a, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        param_ranges(1,j),k10,k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,1));
        title('k_a (absorption)'); hold on;
        xlim([0.7 7])
    end

    if i == 2 % vary k10, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                            k_a,param_ranges(2,j),k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,1)); 
        title('k_1_0 (elimination)');hold on;
        xlim([0.7 7])
    end
    
    if i == 3 % vary k12, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,param_ranges(i,j),k21,V1);
        subplot(2,2,i); plot(t,C(:,1)); 
        title('k_1_2 (distribution)'); hold on;
        xlim([0.7 7])
    end

    if i == 4 % vary k21, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,param_ranges(i,j),V1);
        subplot(2,2,i); plot(t,C(:,1)); 
        title('k_2_1 (redistribution)');hold on;
        xlim([0.7 7])
    end
    
    xlabel('Time (h)'), ylabel('Concentration (ng/M)')
    %suptitle('Sensitivity Analysis for the 4 kinetic rates');
    
    %{
    if i == 5 % vary V1, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,k21,param_ranges(i,j));

        % Plot seperately
        if j==1; figure;end
        plot(t,C(:,1)); hold on;
    end
    %}
    
   end
end

%%%%%%%%%%%%% Plot t,C(:,2) oral dose %%%%%%%%%%%%%
% This could be done in a for loop but it would be hard to deal with 
% all of the formatting so I just broke it up.

figure
for i = 1:5
    for j = 1:7
  
    if i == 1 % vary k_a, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        param_ranges(1,j),k10,k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,2));
        title('k_a (absorption)'); hold on;
        xlim([0.7 10]); ylim([0 7])
    end

    if i == 2 % vary k10, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                            k_a,param_ranges(2,j),k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,2)); 
        title('k_1_0 (elimination)');hold on;
        xlim([0.7 10]); ylim([0 7])
    end
    
    if i == 3 % vary k12, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,param_ranges(i,j),k21,V1);
        subplot(2,2,i); plot(t,C(:,2)); 
        title('k_1_2 (distribution)');hold on;
        xlim([0.7 10]); ylim([0 7])
    end

    if i == 4 % vary k21, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,param_ranges(i,j),V1);
        subplot(2,2,i); plot(t,C(:,2)); 
        title('k_2_1 (redistribution)');hold on;
        xlim([0.7 10]); ylim([0 7])
    end
    
    xlabel('Time (h)'), ylabel('Concentration (ng/M)')
    %title('Sensitivity Analysis for the 4 kinetic rates')
    
    %{
    if i == 5 % vary V1, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,k21,param_ranges(i,j));

        % Plot seperately
        if j==1; figure;end
        plot(t,C(:,2)); hold on;
    end
    %}
    
   end
end


%%%%%%%%% Plot t,C(:,3) the peripheral compartment %%%%%%%%%%%%% 

tspan = 1:0.01:36;

figure
for i = 1:5
    for j = 1:7

    if i == 1 % vary k_a, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        param_ranges(1,j),k10,k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,3)); 
        title('k_a (absorption)');hold on;
        xlim([0.7 21]); ylim([0 8]) 
    end

    if i == 2 % vary k10, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                            k_a,param_ranges(2,j),k12,k21,V1);
        subplot(2,2,i); plot(t,C(:,3));
        title('k_1_0 (elimination)');hold on;
        xlim([0.7 35]); ylim([0 8])
    end
    
    if i == 3 % vary k12, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,param_ranges(i,j),k21,V1);
        subplot(2,2,i); plot(t,C(:,3)); 
        title('k_1_2 (distribution)'); hold on;
        xlim([0.7 35]); ylim([0 8])
    end

    if i == 4 % vary k21, hold other params constant
        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,param_ranges(i,j),V1);
        subplot(2,2,i); plot(t,C(:,3)); 
        title('k_2_1 (redistribution)');hold on;
        xlim([0.7 35]); ylim([0 8])
    end
    
    xlabel('Time (h)'), ylabel('Concentration (ng/M)')
    %title('Sensitivity Analysis for the 4 kinetic rates')
    
    %{
    if i == 5 % vary V1, hold other params constant

        [t, C] = ode45('C_single_dose', tspan, C0,[],SA,...
                        k_a,k10,k12,k21,param_ranges(i,j));

        % Plot seperately
        if j==1; figure;end
        plot(t,C(:,2)); hold on;
    end
    %}
    
   end
end
