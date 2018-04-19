
% Ben Stear
% Quatitative Systems Biology
% Yeast Micro Array Analysis
% 4/4/18

%% Download and Parse GSE and GPL data

gse = geoseriesread('/Users/dawnstear/downloads/GSE17257_series_matrix.txt');
gpl = bmes_downloadandparsegpl('GPL2529');       
get(gse.Data);
d = gse.Data;
 
gpl_genes = gpl.Data(:,11);   % 'gene symbol' 
gpl_probes = gpl.Data(:,1);   % gpl probe ID
gse_probes = d.rownames;       % gse probe ID

map = containers.Map(gpl_probes,1:numel(gpl_probes));
for i=1:numel(gse_probes)
	if map.isKey(gse_probes{i}); 
        MAP_GSE_GPL(i)=map(gse_probes{i}); end
end

gse_genes = gse_probes; %make a copy, so entries not found will keep the probe name.
gse_genes(find(MAP_GSE_GPL)) = gpl_genes(MAP_GSE_GPL(find(MAP_GSE_GPL)));

table(gse_probes(1:5), gse_genes(1:5),'VariableNames',{'gseprobe','gsegene'});
d = d.rownames(':',gse_genes); % now replace gseprobes with gsegenes

d_array = single(d); % convert to numerical array 


%% Data Analysis: Find differentially expressed genes between groups of samples.

% The data in the datamatrix d ranged from ~0 - ~1000 which would indicate
% that the data has not been log2 normalized, so I did that first.
d_array = log2(d_array); 

% Then I grabbed the first 3 columns which correspond to the DSMO
% group and the last 3 columns which correspond to the CQ group
DMSO = d_array(:,1:3);   
CQ   = d_array(:,4:6); 

% From the paper "The fold changes were calculated by ratio 
% of signals in CQ-treated samples to that in DMSO-treated 
% controls, and presented as the averages of three experiments."
DMSOave = sum(DMSO,2)/3;
CQave = sum(CQ,2)/3;

% subtract average expression from the two groups. If they were not
% log2 normalized we'd have to divide the two groups.
log2fc = CQave - DMSOave;

% k = [];
% g = 3.77;
% while (numel(k)~= 11)      % find the top 11 genes with the highest 
% k = find(abs(log2fc) > g); % abs(fc). I did 11 and not 10 because
% g =- 0.001;                % one of the genes I found in the top 10
% end                        % had a blank geneID


k = find(abs(log2fc) > 3.760); % use the cutoff value I found, 3.760
topfc = zeros(1,numel(k));     % and find the top 11 genes directly
                               % instead of using the loop every time
for i = 1:numel(k)
    topfc = d.rownames{k(i)};   % print top 10 (11) genes
end

for i = 1:numel(k)
    log2fc(k);
end
% find pvals, and then correct for false positives
[~,pvals]=ttest2(DMSO', CQ');
fpvals = mafdr(pvals);
[dpvals]=mattest(DMSO,CQ , 'permute',10);

%% Calculate fold change
% convert log2fc to negfc
negfc = 2.^log2fc;
negfc(negfc<1) = - 1./negfc(negfc<1);
negfc_sorted = sort(negfc);
dpvals=[dpvals bioma.data.DataMatrix(negfc,'ColNames',{'negfc'})];


%% Data Analysis: Gene Set Enrichment

gplgobio  = gpl.Data(:, strcmp(gpl.ColumnNames,...
                     'Gene Ontology Biological Process'));

% Clean up Gene Ontology terms .
for i=1:numel(gplgobio)
	if isempty(gplgobio{i}); gplgobio{i}={};
	else gplgobio{i} = strsplit(gplgobio{i},' /// '); end
end

GOBIO = unique( [gplgobio{:}] );
map = containers.Map(GOBIO,1:numel(GOBIO));
for i=1:numel(gplgobio)
	for j=1:numel(gplgobio{i})
		gplgobio{i}{j} = map(gplgobio{i}{j});
	end
	gplgobio{i} = cell2mat(gplgobio{i});
end

gsegobio = cell(size(d,1),1); 
gsegobio(find(MAP_GSE_GPL)) = gplgobio(MAP_GSE_GPL(find(MAP_GSE_GPL)));

% Find the terms that have at least one significant gene for it.
Isig=dpvals(:,'')<=0.01;
dsig = d(Isig,:);
dpvalssig = dpvals(Isig,:);
gsegobiosig = gsegobio(Isig);
candidategobio = unique( [gsegobiosig{:}] );
gobioid = candidategobio(1);
GOBIO{ gobioid };

x = nnz([gsegobiosig{:}] == gobioid);
K = nnz([gsegobio{:}] == gobioid);
N = nnz( ~cellfun(@isempty,gsegobiosig) );
M = nnz( ~cellfun(@isempty,gsegobio) );

pval = hygecdf(x,M,K,N,'upper'); %pval = 1-hygecdf(x,M,K,N);

% Gene Set Enrichment is reported in the word document.
