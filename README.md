SDpop
=====

Author: Jos Kafer, jos.kafer@cnrs.fr

SDpop infers sex-linkage from genotyping data of several individuals of both sexes, collected in panmictic populations. 

To install, make sure you have a sufficiently recent C++ compiler (understanding C++11) and type:
make all
If you're on a computer cluster, you might want to compile statically, which you could try to do with:
make CXXFLAGS="--static" all

# Preamble

The software is provided "as is" and is still under development. 
It is developed by me (Jos Käfer); I'm a research scientist interested in plant evolutionary biology, not a software 
developer interested in optimizing code. Thus, the programs provided here are aimed to work on my data and with my 
pipelines. I've tried to describe the data formats and caveats, but certainly not in an exhaustive manner. Error 
handling in the code is rudimentary, so running the programs on you data will very likely cause segmentation faults (or 
worse).

You are free to use the code, as long as you cite the paper in your work. If you run into any issue, and you're not 
familiar enough with C/C++ to solve it yourself (or not enduring enough to scroll through my thousands of lines of 
ill-commented code, which I understand), feel free to contact me. I'm happy to make the code more robust, and learn 
something about your research. 

If you're a software developer interested in optimizing the code, I welcome your help and suggestions. The code can be 
written more concisely, and several algorithmic techniques could greatly diminish the running time and memory usage. One
of the first things to do would be to compress the data and to do calculations for sites with the same observed 
genotypes only once.

# Typical data analysis workflow

You can read the paper about SDpop here: https://academic.oup.com/genetics/article/218/2/iyab053/6187941
It presents the analyses of two datasets.

I'll summarize the steps used in a typical analysis:

0) Get samples from a species, preferably from a single population. If you don't know what to expect, 10 female and 10 male individuals might be a good start.
Sequence them using an appropriate technique. There are many techniques available, and SDpop can work with RNAseq data, DNAseq, RADseq, and probably even more.

1) Map the sequence reads on a reference to obtain a vcf file covering all positions (sometimes called a gvcf file). 
There are many tools available, depending on your data, reference, etc. This step is crucial and usually takes most time and
energy in the analysis. There are some considerations in the paper and there is some practical help in the user guide.

2) Check for unwanted population structure. I recommend doing a PCA on genetic variation, for example with plink (https://www.cog-genomics.org/plink2/) and 
discarding individuals that are clearly not interbreeding with the others. You could also check whether your individuals are siblings or clones, for example using 
"relatedness2" from vcftools. All biologists know that panmixia is an ideal state never reached, but the closer your sample is to panmixia, the more reliable the results will be.

3) If needed, use SDpop to test whether the species actually has sex chromosomes. Run the models without sex chromosomes, with XY chromosomes, and with ZW chromosomes. 
The "XY and ZW" model is provided for special cases and is not biologically relevant in most cases, so I recommend not spending too much time on figuring out why 
it has a lower BIC than the XY or the ZW models. 

4) Once you've figured out whether your species has XY or ZW chromosomes, identify the sex-linked genes or regions. If you're working on genes or small contigs, a score of 0.8 
would be a first threshold to explore. Note that the detection of XY (or ZW) genes or contigs is much more robust than the detection of X- (or Z-) hemizygous genes or contigs. 
Read the paper to understand why.

5) Predict the X and Y (or Z and W) gametolog sequences and study their divergence.

6) Carefully check the inferences using as many other pieces of information you have. If you have a well-assembled genome to place your contigs or genes on, you're lucky. 
For the X- or Z-hemizygous genes or regions, coverage information is almost certainly needed to confirm the inferences.



