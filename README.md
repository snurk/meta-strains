# meta-strains


## To install modified Clomial package

In R:

* install.packages("matrixStats", lib="\~/R/libs")

* install.packages("permute", lib="\~/R/libs", repos='http://cran.rstudio.com/')

In console:

* R CMD INSTALL -l ~/R/libs/ Clomial_meta


## To run Clomial pipeline

1. Create an empty directory in "experiments".

2. Copy your reference file to directory "refs" and run "bwa index" on it.

3. Copy "Snakefile_Clomial" to your directory and rename it to "Snakefile".

4. Change settings in first lines of your Snakefile.

5. Run "snakemake --cores (number_of_threads)".

6. Do not forget to clear GitHub repository from gigabytes of data.


## To run PhyloWGS pipeline

Same as Clomial pipeline, but for "Snakefile_PhyloWGS"
