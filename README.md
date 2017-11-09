# MetagenomicDC

### MetagenomicDC: Metagenomic Data Classifier based on deep learning models
A 16S short-read sequences classification tool based on k-mer representation and deep learning architecture (CNN and DBN).
It can be integrated into the most common pipelines for metagenomic analysis.
It can be successfully used for classifying both SG and AMP data.

#### Dependences:
* Python (2.7.x)
* Theano (0.8.2)
* Keras library (2.x)


###
### Datasets
This distribution contains three datasets. They are available in "__data__" folder. 
* The "__1000seq.fasta__" dataset has been used for experiments in the manuscript submitted for pubblication at BMC Bioinformatics journal. It is composed by a training dataset with 1000 16S fasta sequences in unaligned fasta format, with 100 genera and 10 species of each genus. All of them are extracted from Rfam database (release 11, update 5 dated September 30, 2016). 
* The "__16S-AMP-trimmed.fa.zip__" dataset contains simulated amplicon technology short-reads (with the Grinder tool in V3-V4 hypervariable region using the following primers: ``CCTACGGGAGGCAGCAG`` and ``CCGTCAATTCMTTTRAGT``), from "__1000seq.fasta__" dataset.
* The "__16S-SG-reads.fa.zip__" dataset contains simulated shotgun short-reads (with the Grinder tool), from "__1000seq.fasta__" dataset.

