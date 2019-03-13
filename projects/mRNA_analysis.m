%  Mohammad Bilash Hossain & Ben Stear           %
%                                                %
% Quantitative Systems Biology Final Project     %
%                                                %
%                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Project Title: MicroRNA-135b regulates ER, AR and HIF1AN and 
               affects breast and prostate cancer cell growth
%}


% load data
%  There are 2 groups, the pre miR135b group and the negative control,
%  Scramble group, both groups have 3 time points (12h,24h,36h) and each
%  time point has 2 samples (b1,b2)

gse = bmes_downloadandparsegse('GSE57820');
gpl = bmes_downloadandparsegpl('GPL10558');  
d = gse.Data;

% Pre-miRNA data is in columns 1-6, find average of the b1 and b2 samples
% by adding and dividing by two.
miR135b_12h = (double(d(:,1)) + double(d(:,2)))/2;
miR135b_24h = (double(d(:,3)) + double(d(:,4)))/2;
miR135b_36h = (double(d(:,5)) + double(d(:,6)))/2;

% Same thing for the Scramble negative control data
Scr_12h = (double(d(:,7))+double(d(:,8)))/2; 
Scr_24h = (double(d(:,9))+double(d(:,10)))/2;
Scr_36h = (double(d(:,11))+double(d(:,12)))/2;

% find pvals, and then correct for false positives
x = double(d(:,1:6)); y = double(d(:,7:12));
[~,p]=ttest2(x', y'); p=p';    % pval< 0.05 was used
[FDR, Q] = mafdr(p); Q=Q'; % Qval< 0.05 was used

% Find fold change, then do log2(miR135b) 
fc12 = log2(miR135b_12h ./ Scr_12h) ;% 12 Hours
fc24 = log2(miR135b_24h ./ Scr_24h);% 24 Hours
fc36 = log2(miR135b_36h ./ Scr_36h);% 36 Hours

% Replace any fold change value that is less than 1, with its negative inverse (-1/x).
for i = 1:length(fc12)
    if fc12(i)<1; fc12(i) = -1/fc12(i); end
end
for i = 1:length(fc24)
    if fc24(i)<1; fc24(i) = -1/fc24(i); end
end
for i = 1:length(fc36)
    if fc36(i)<1; fc36(i) = -1/fc36(i); end
end

% Find indices where p& > 0.05, set that value in the fc = 0 bc we dont want it
fc12([p>0.0500]) = 0; 
fc24([p>0.0500]) = 0; 
fc36([p>0.0500]) = 0; 

% Find expression values that greater or less than |5000|
fc12_up   = fc12 >  5000;   fc24_up   = fc24 >  5000;   fc36_up   = fc36 >  5000; 
fc12_down = fc12 < -5000;   fc24_down = fc24 < -5000;   fc36_down = fc36 < -5000;

% % Find which the gene probes that correspond to the differentially
% expressed indices.
geneprobes = (d.RowNames);
genenames = gpl.Data(:,13);

fc12_up_probes = geneprobes(fc12_up);
fc12_down_probes = geneprobes(fc12_down);
fc24_up_probes = geneprobes(fc24_up);
fc24_down_probes = geneprobes(fc24_down);
fc36_up_probes = geneprobes(fc36_up);
fc36_down_probes = geneprobes(fc36_down);

m=1;
for i = 1:length(fc12_up_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc12_up_probes(i)) == cell2mat(geneprobes(j))
            fc12up_genes(m) = string(genenames(j));
            m=m+1;
        end
    end
end
% delete the 6th entry, no gene ID was found 
fc12up_genes(6) = [];

m=1;
for i = 1:length(fc12_down_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc12_down_probes(i)) == cell2mat(geneprobes(j))
            fc12down_genes(m) = string(genenames(j)); m=m+1;
        end
    end
end
% delete the 1st entry, gene ID was 'permuted negative'
fc12down_genes(1) = [];

m=1;
for i = 1:length(fc24_down_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc24_down_probes(i)) == cell2mat(geneprobes(j))
            fc24down_genes(m) = string(genenames(j)); m=m+1;
        end
    end
end

m=1;
for i = 1:length(fc24_up_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc24_up_probes(i)) == cell2mat(geneprobes(j))
            fc24up_genes(m) = string(genenames(j)); m=m+1;
        end
    end
end

m=1;
for i = 1:length(fc36_up_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc36_up_probes(i)) == cell2mat(geneprobes(j))
            fc36up_genes(m) = string(genenames(j)); m=m+1;
        end
    end
end

m=1;
for i = 1:length(fc36_down_probes)
    for j = 1:length(geneprobes)
        if cell2mat(fc36_down_probes(i)) == cell2mat(geneprobes(j))
            fc36down_genes(m) = string(genenames(j)); m=m+1;
        end
    end
end

% Find genes that upregulated by all three groups
all_up = sum(fc12_up&fc24_up&fc36_up);  % no genes are in all 3 upregulated time samples
all_down = sum(fc12_down&fc24_down&fc36_down); % also didnt find any genes across the 3 groups

% Genes that are upregulated by >5000x
fc12up_genes = ['LOC650518','C1orf85', 'LOC645743','ARHGAP6','MEF2D','KRT79','UBE2A','RAB27A','TATDN3'];
fc24up_genes = ['MDH1', 'AP4S1','PDE4A','LOC728002'];
fc36up_genes = ['LOC728116','LOC100129928'];

% Genes that are downregulated by >5000
fc12down_genes = [ 'SEMA3B', 'BOLA2', 'CDH8','DTWD2','SSX4'];
fc24down_genes = ['LOC728393','SCARNA10','GHRLOS'];
fc36down_genes = ['LOC643995', 'TARP','CXorf23','IPP'];

% none of these^^^ genes matched any of the genes in fig 5c from the paper


% Make Venn Diagrams showing the overlap of genes differentially
% expressed at 12h, 24h and 36h. Use a smaller threshold than 5000 to find
% differentially expressed genes.
fc12_up   = fc12 >  50;   fc24_up   = fc24 >  50;   fc36_up   = fc36 >  50; 
fc12_down = fc12 < -50;   fc24_down = fc24 < -50;   fc36_down = fc36 < -50;

% Find all over or under expressed genes for each time sample
fc12_tot = (fc12_up | fc12_down); % 3714 genes
fc24_tot = (fc24_up | fc24_down); % 3566 genes
fc36_tot = (fc36_up | fc36_down); % 3496 genes

% % For a 3-circle venn, Z is a 7 element vector [z1 z2 z3 z12 z13 z23 z123]
% % Compute unions for Venn Diagram
z1 = sum(fc12_tot);  
z2 = sum(fc24_tot);  
z3 = sum(fc36_tot);  
z12 = sum(fc12_tot&fc24_tot); z23 = sum(fc24_tot&fc36_tot);
z13 = sum(fc12_tot&fc36_tot); z123 = sum(fc12_tot&fc24_tot&fc36_tot);

% % Create Venn Diagram
venn([z1 z2 z3 z12 z23 z13 z123])

% red circle = 12h, green circle = 24h, blue circle = 36h
% 43 genes are up or down regulated by 50x across all 3 time samples

% Find indices of these 43 genes
idx = find(fc12_tot&fc24_tot&fc36_tot);
% Find gene probes for these 43 genes
genes43_probes = geneprobes(idx);
m=1; 
for i = 1:length(genes43_probes)
    for j = 1:length(geneprobes)
        if cell2mat(genes43_probes(i)) == cell2mat(geneprobes(j))
            genes43(m) = string(genenames(j));
            m=m+1;
        end
    end
end

disp(genes43) % This list has several in common with the 64 
              % genes the paper listed that they found up/down regulated
              % at all 3 time samples



