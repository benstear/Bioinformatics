

% plasma concentration times over 4 wks
% Data Source: 'Paroxetine serum concentrations in 
%               depressed patients and response to treatment.'


Cplasma = [ 19 20 22.7 32.2 40 43 45 48.8 53.4 45 39.1 36 34];
t = 1:length(Cplasma);
sz = 55;
% scatter(x,Cplasma,sz,'filled')



% insert our data and convert to table
data = array2table([t; Cplasma]');

%data = array2table(data);
data.Properties.VariableNames = {'Time' 'Conc'};

% put in grouped data format for smbiofit()
gData = groupedData(data);
gData.Properties.VariableUnits = {'hour','milligram/milliliter'};

pkm = PKModelDesign;
pkc1 = addCompartment(pkm, 'Central', 'DosingType', 'Bolus', ...
                     'EliminationType', 'Linear-Clearance'); % enzymatic
                 
model = construct(pkm);

configset = getconfigset(model);
%configset.CompileOptions.UnitConversion = true;

dose = sbiodose('dose');
dose.TargetName = 'Drug_Central';
dose.StartTime = 0;
dose.Amount = 40;
dose.AmountUnits = 'milligram';
dose.TimeUnits = 'hour';

responseMap = {'Drug_Central = Conc'};
paramsToEstimate = {'Central','Cl_Central'};
estimatedParams = estimatedInfo(paramsToEstimate,'InitialValue',[1 1],'Bounds',[1 5;0.5 2]);

fitConst = sbiofit(model,gData,responseMap,estimatedParams,dose);
%{
fitProp = sbiofit(model,gData,responseMap,estimatedParams,dose,...
                      'ErrorModel','proportional');
fitExp  = sbiofit(model,gData,responseMap,estimatedParams,dose,...
                      'ErrorModel','exponential');
%}
s.Labels.XLabel = 'Time (hour)';
s.Labels.YLabel = 'Concentration (milligram/milliliter)';
plot(fitConst,'AxesStyle',s);




