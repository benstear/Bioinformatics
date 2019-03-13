
D = 30;
F = 0.5;
k_a = 14;
k10 = 2.9;
Vd = 8.7;
tspan = 1:301;
C = zeros(1,length(tspan));
c(1) = 0;

for i = 2:length(tspan)
    C(i) = ((F*D*k_a)/Vd*(k_a-k10))*(exp(-k10*tspan(i)))-(exp(-k_a*tspan(i)));
end