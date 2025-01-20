# LVgs: A robust K-mer-based genome size evaluation pipeline
LVgs is a novel strategy for estimating the **g**enome **s**ize robustly by plotting iteration of the large-scale estimates through continuous‌ ascent of K-mer length and searching the **L**imiting **V**alue of genome forecast with taking advantage of HiFi reads.

## Installation
LVgs is shell-based and R-based scripts, which is compile-free.

### 1. Pre-requisites
The following requirements are assumed to be installed (with path to their executables available in $PATH).

- [FastK](https://github.com/thegenemyers/FASTK)
- [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)
- [HiFiasm](https://github.com/chhylp123/hifiasm)
- [tseries](https://CRAN.R-project.org/package=tseries)

### 2. Download
```
git clone https://github.com/xingjianfeng100/LVgs.git
cd LVgs
```
or
```
wget https://github.com/xingjianfeng100/LVgs/releases/download/*/LVgs-***.tar.gz
tar -xzvf  LVgs-***.tar.gz
cd LVgs-*** 
```
### 3. Set path to LVgs available in environment
```
export PATH=$PWD:$PATH
```
## Usage

```console
Usage: 	LVgs.sh -b stat_k-mer_size -e end_k-mer_size -i step_size -r reads
-b <INT>                                      K-mer length for beginning of loop forecasts (Bp)                                      [required]
-e <INT>                                      K-mer length for ending of loop forecasts (Bp)                                         [required]
-s <INT>                                      K-mer length iteration step size for loop forecasts (Bp)                               [required]
-r  <.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz]> datasets for generating K-mers                                                         [required] 
-p  <INT>                                     ploidy (1, 2, 3, 4, 5, or 6) for GenomeScope2 to use                                   [default: 2]
-f                                            full mode to keep intermediate files                                                   [default: disabled]
-c                                            re-correct the HiFi consensus reads using hifiasm                                      [default: disabled]
-o  <output_prefix>                           prefix of the output files                                                             [default: out]
-t  <INT>                                     the number of threads to use                                                           [default: 1]
-M  <INT>                                     Use -M GB of memory in K-mer counting steps of FastK                                   [default: 12]
-W  <workpdir>                                directory to execute loop forecasts                                                    [default: none]
-k <"-v|-t|-p|-bc|-c">                        other options for FastK, details see https://github.com/thegenemyers/FASTK             [default: none]
-g <"-l|-m......">                            other options for GenomeScope2, details see https://github.com/tbenavi1/genomescope2.0 [default: none]
-V                                            show version number
-h                                            display this help and exit
```
> [!NOTE]
> * ```-r <.cram|.[bs]am|.db|.dam|.f[ast][aq][.gz]>``` inputting reads need to specify absolute path.
> * Other options for FastK or GenomeScope2 should be enclosed in quotes, e.g., ```-k "-v -t4" -g "-m 1000000 --verbose"```.

### Outputs

```console
*.report				#Final evaluation report for LVgs
*.stat.tab				#Forecasts statistics with multi-K
*figure/				#Genome size (GS) curve as well as other seven auxiliary indicators, Forecasted by multi-K
    ├── GS_linear_plot.pdf
    ├── kcov_plot.pdf
    ├── kmer_total_plot.pdf
    ├── kmer_type_plot.pdf
    ├── Model_fit_plot.pdf
    ├── Read_Error_plot.pdf
    ├── Repeat_linear_plot.pdf
    |__ Unique_linear_plot.pdf
*_limiV/				#K-mer spectrum and K-mer histogram files with final-evaluation K
    ├── kmer.hist			#K-mer histogram file produced by FastK
    ├── k*.his				#A special ASCII K-mer histogram file suitable for input to GenomeScope2 produced by Histex from 'kmer.hist' 
    |$\underline{  }$__ gsco_k*/			#K-mer spectrums generated by GenomeScope2
```


## Examples
### Demo data downloading
A HiFi dateset generated from _Phthorimaea absoluta_ was used as an example.
```
wget  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR234/030/SRR23497930/SRR23497930_subreads.fastq.gz
```
### Quick Start
```
LVgs.sh -b 17 -e 77 -s 10 -r /absolute/path/SRR23497930_subreads.fastq.gz  -t 32 -W demo
```
result files are stored in ```demo/```


## Contact
Other than raising issues in github, you can contact xingjianfeng100@126.com for specific issues.

## Cite
If you use our tools, please also cite these remarkable tools we have invoked in LVgs:
+    FastK: https://github.com/thegenemyers/FASTK
+    GenomeScope2: Ranallo-Benavidez, T. R., Jaron, K. S. & Schatz, M. C. GenomeScope 2.0 and Smudgeplot for reference-free profiling of polyploid genomes. Nat. Commun. 11, doi:10.1038/s41467-020-14998-3 (2020)
+    HiFiasm: Cheng, H., Asri, M., Lucas, J., Koren, S. & Li, H. Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. Nat Methods 21, 967-970, doi:10.1038/s41592-024-02269-8 (2024).
+    tseries: Trapletti, A. & Hornik, K. tseries: Time Series Analysis and Computational Finance. (R package version 0.10-58: https://CRAN.R-project.org/package=tseries, 2024)
