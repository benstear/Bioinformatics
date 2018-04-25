% Ben Stear
%
%  MicroRNA analysis
%
% Analyze miRNA levels in patient w & w/o Complex Regional Pain Syndrome (CRPS)

% Use the unfiltered Excel data.
unfilt = xlsread('/Users/dawnstear/desktop/CRPS_unfiltered.xlsx');
numrows = size(unfilt,1);
numcols = size(unfilt,2);
rows2delete = zeros(1,150);
h = 1;

% find rows with fewer than 3 numerical entries (60 or more NaNs) so we can delete them
for i = 1:numrows
      nancount = find(isnan(unfilt(i,:))); % count how many NaNs per row
      if length(nancount)>58  % length(nancount) == 60 || 61 means only 1 or 2 non-NaN terms
          rows2delete(h) = i; % save rows we want to delete
          h = h+1;
      end
 end
          
 % Now delete those rows  
 unfilt2 = unfilt; % save original data
 unfilt2([rows2delete],:) = []; % delete those rows
 
numrows = size(unfilt2,1);
rowaves = zeros(numrows,1);
rowsums = nansum(unfilt2([1:581],:),2); % use nansum to find sums of rows while ignoring NaNs
numcount = 0; 

for i = 1:numrows 
   numcount = sum(~isnan(unfilt2(i,:))); % find out how many values we have (ie Non-NaN terms)
   rowaves(i) = rowsums(i)/numcount; % compute row average
   numcount = 0; % reset numcount for next row
end

% if index == NaN replace that whole with ith rowaverage
for i =1:numrows
    for j = 1:numcols
        if isnan(unfilt2(i,j))
            unfilt2(i,j)= rowaves(i);
        end
    end
end

%%------------------------------------------------------------------------%
% Now Use the filtered Excel data.
filt = xlsread('/Users/dawnstear/desktop/CRPS_filtered.xlsx');

% Find column average for each control sample
MammU6  = filt(87,:); 
RNU48   = filt(96,:); 
RNU48_2 = filt(225,:);
RNU44   = filt(184,:); 
RNU44_2 = filt(280,:);
controls = [88,97,226,185,281] - 1; % subtract by 1 b/c row #'s are shifted up 1 in matlab matrix compared to .xlsx file 
colaves = [1,size(filt,2)];         

for i=1:size(filt,2) % #cols = (61)
    colaves(i) = sum((filt(controls,i)))/numel(controls); % calculate average of controls for each col
end
                                                                                                                                        
 % normalized deltaCT values
filtnorm = filt-colaves; 


% low ct value = highly expressed
% ct stands for cycle threshold and represents how many cycle of pcr we
% need to do to get a detectable lvl of that miRNA

% mean(A,2) = average by each row, here, rows = miRNAA; so we find ave. of
% each miRNA expression level. these are average deltaCT's for each respective group
filtControls  = mean(filt(:,1:20),2); 
filtCRPS      = mean(filt(:,21:61),2);

% deltadeltaCT = filtControls - filtCRPS
ddCT = filtControls - filtCRPS;  

% Calculate 2^-deltadeltaCT, which are the fold change values.
ddCT_fc = 2.^(-ddCT); 

% Replace any fold change value x that is less than 1, with its negative inverse (-1/x). 
for i = 1:numel(ddCT_fc)
    if ddCT_fc(i) < 1
        ddCT_fc(i) = -1/ddCT_fc(i);
    end
end

index = find(ddCT_fc> 2.2755); % get index of top10 highest fold changes
toptenfc = ddCT_fc(index);     % and use them to get the actual fc vals


% Find p-values
% fpvals = mafdr(pvals);
[~,pvals] = ttest2(filt(:,1:20)',filt(:,21:61)'); %%%%% pvals on filt or ddCT_fc ???????
pvals = pvals';
Fold_Change= [2.4067;2.4387; 2.8309; 2.2761; 3.0093; 2.3874; 2.4953; 2.2755; 3.9236; 4.4489];
miRNAs = {'hsa-miR-218';'hsa-miR-100'; 'hsa-let-7e';'hsa-miR-361-5p';...
          'hsa-miR-574-3p';'hsa-miR-31';'hsa-miR-192';'hsa-miR-18a#';...
          'hsa-miR-130b#';'hsa-miR-664'};
p_value = [22.1963e-003; 46.0825e-003; 18.7397e-003; 7.0276e-003;...
           293.3798e-003; 72.6989e-006; 20.1732e-006; 80.6501e-003; 25.0621e-006; 37.6428e-006];

% create table with miRNA name, fold changes and p-vals       
microRNA = table(miRNAs,Fold_Change,p_value);

% sort by p-value (3rd column)
microRNA = sortrows(microRNA,3);

% Now use targetscan.org to find mRNA targets of these microRNA. 
% Then use DAVID to perform gene enrichment and functional annotation




