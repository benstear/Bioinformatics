

% 2 compartment model with sim bio

load(fullfile(matlabroot,'examples','simbio','data10_32R.mat'))

gData = groupedData(data);
gData.Properties.VariableUnits = {'','hour','milligram/liter','milligram/liter'};

sbiotrellis(gData,'ID','Time',{'CentralConc','PeripheralConc'},...
            'Marker','+','LineStyle','none');
        
pkmd                 = PKModelDesign;
pkc1                 = addCompartment(pkmd,'Central');
pkc1.DosingType      = 'Infusion';
pkc1.EliminationType = 'linear-clearance';
pkc1.HasResponseVariable = true;
pkc2                 = addCompartment(pkmd,'Peripheral');
model                = construct(pkmd);
configset            = getconfigset(model);
configset.CompileOptions.UnitConversion = true;


dose             = sbiodose('dose','TargetName','Drug_Central');
dose.StartTime   = 0;
dose.Amount      = 100;
dose.Rate        = 50;
dose.AmountUnits = 'milligram';
dose.TimeUnits   = 'hour';
dose.RateUnits   = 'milligram/hour';

responseMap = {'Drug_Central = CentralConc','Drug_Peripheral = PeripheralConc'};

paramsToEstimate   = {'log(Central)','log(Peripheral)','Q12','Cl_Central'};
estimatedParam     = estimatedInfo(paramsToEstimate,'InitialValue',[1 1 1 1]);