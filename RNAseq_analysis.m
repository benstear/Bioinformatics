
% Ben Stear                             %
%                                       %
% RNA sequence analysis                 %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
In this analysis we looked at RNA sequencing data. The data was read from the sequencer into fastq file format. 
We used Burrows Wheeler Alignment to first index the yeast reference genome using the command $ bwa index ref.fa .

Paper the code is based on: Regulation of Glucose-Dependent Gene Expression by the RNA Helicase Dbp2 in Saccharomyces cerevisiae

yeast type: saccharomyces cerevisiae 
genbank ncbi BY4741

dbp2-mutant: SRR1302792
wild type:   SRR1302790

File pipeline:
.fastq --> .SAM --> .txt --> .m
%}


%% BWA

%{
What I entered into cmd line:

git clone https://github.com/lh3/bwa.git
cd bwa; make
./bwa index ./genome.fa
./bwa mem genome.fa ./SRR1302792_pass.randsample.fastq > mapreads92.sam
./bwa mem genome.fa ./SRR1302792_pass.randsample.fastq > mapreads90.sam

%}


%% featureCounts 

%{
featureCounts -a <annotation_file> -o <output_file> input_file1

What I entered into the cmd line:

For dbp2-MUTANT:
./featureCounts -a /Users/dawnstear/Desktop/RNA_QSB/genes.gtf 
        -o counts92.txt /Users/dawnstear/bwa/mapreads92.sam

For WILD TYPE:
./featureCounts -a /Users/dawnstear/Desktop/RNA_QSB/genes.gtf 
        -o counts90.txt /Users/dawnstear/bwa/mapreads90.sam
%}

% Now read in as a table
ftrCnts90 = readtable('/Users/dawnstear/Desktop/RNA_QSB/counts90.txt');
ftrCnts92 = readtable('/Users/dawnstear/Desktop/RNA_QSB/counts92.txt');

% Isolate gene expression data (in 7th col), convert table2array, add
% 1 so we dont get division by 0 error when calculating fold change
geneEx90 = table2array(ftrCnts90(:,7))+1;
geneEx92 = table2array(ftrCnts92(:,7))+1;

% Divide the dbp2-mutant value (SRR1302792) by the wild type
% value (SRR1302790) to find fold change.
foldchange = geneEx92(:,1)./geneEx90(:,1);

% Need to make values less than 1 negative, we can do: -1/n  
for i =1:length(foldchange)
   if foldchange(i)<1
       foldchange(i) = (-1 / foldchange(i));
   end
end

% Shouldve done something like this: 
%    foldchange = -1./ foldchange(foldchange<1);


% Because the top10 absval fold changes are all >111, I used 111 as
% a threshold to find the indices I need in order to find the names
% of the genes.
idx = find(abs(foldchange) > 111);

% The genes are in the same order in both tables, so use either one
Gene = table2array(ftrCnts90(idx,1));
Fold_Change = foldchange(idx);

% Create table 
T = table(Gene,Fold_Change)

%% Discussion

%{

Because I divided mutant by wildtype, a positive fold change in a gene 
would mean it was expressed more highly in the mutant, likewise a 
negative fold change means a gene was expressed less in the mutant
compared to wild type.

The only gene that made the top ten in absolute value fold change that 
was negatively expressed was 'YNL112W' which coded for the DEAD-Box 
protein dbp2. This makes sense because this is the protein that the
paper mainly focused on. My results that dbp2 is more highly expressed
in the wild type are in agreement with the paper, "we report that
either a carbon source switch or glucose deprivation results in rapid 
export of Dbp2 to the cytoplasm."

By far the gene that had the highest fold change was 'YFL014W' (1,463x 
higher in mutant) which  codes for the heat shock protein HSP12. 
HSP12 is a plasma membrane protein which is involved in "maintaining 
organization during stress conditions; induced by heat shock, oxidative 
stress, glucose depletion, oleate and alcohol; protein abundance 
increased in response to DNA replication stress and dietary restriction; 
regulated by the HOG and Ras-Pka pathways; required for dietary 
restriction-induced lifespan extension." This description is also in 
agreement with the paper and its findings.

%}

