`rlsim` - a package for simulating RNA-seq library preparation with parameter estimation
========================================================================================

![rlsim architecture](https://raw.github.com/sbotond/rlsim/master/doc/pix/rlsim_main.png)

What is the `rlsim` package?
--------------------------

The `rlsim` package is a collection of tools for simulating [RNA-seq](http://en.wikipedia.org/wiki/RNA-Seq)
library construction, aiming to reproduce the most important factors
which are known to introduce significant biases in the currently used
protocols: [hexamer priming](http://dx.doi.org/10.1093/nar/gkq224), 
[PCR amplification](http://dx.doi.org/10.1093/nar/gks001) and [size selection](http://dx.doi.org/10.1093/nar/gkq1015).
It allows for a systematic exploration of the effects of the individual biasing
factors and their interactions on downstream applications by simulating
data under a variety of parameter sets.

The implicit simulation model implemented in the main tool (`rlsim`) is
inspired by the actual library preparation protocols and it is more
general than the models used by the bias correction methods hence it allows for 
a fair assessment of their performance.

Although the simulation model was kept as simple as possible in order to
aid usability, it still has too many parameters to be inferred from data
produced by standard RNA-seq experiments. However, simulating datasets
with properties similar to specific datasets is often useful. To address
this, the package provides a tool (`effest`) implementing simple
approaches for estimating the parameters which can be recovered from
standard RNA-seq data (GC-dependent amplification efficiencies, fragment
size distribution, relative expression levels).

[![Catalogued on GSR](http://popmodels.cancercontrol.cancer.gov/static/img/gsr_tile.jpg)](http://popmodels.cancercontrol.cancer.gov/gsr/packages/rlsim)

Citing the `rlsim` package
--------------------------

An associated manuscript is in preparation, meanwhile the package should
be cited as:

 - Botond Sipos, Greg Slodkowicz, Tim Massingham, Nick Goldman (2013) *Realistic simulations reveal extensive sample-specificity of RNA-seq biases* *arXiv*:[1308.3172](http://bit.ly/rlsax) 

The analysis pipeline used to generate the results is available at [github.com/sbotond/paper-rlsim](http://bit.ly/rlsim-pl).

`rlsim` was brought to you by the [Goldman group](http://www.ebi.ac.uk/research/goldman) from [EMBL-EBI](http://www.ebi.ac.uk).

Key features
------------

-   Simulation of priming biases loosely based on a nearest-neighbor
    thermodynamic model.

-   Exact simulation of PCR amplification on the level of individual
    fragments (consistent across expression levels, no approximations).

-   Fragment-specific amplification efficiencies determined by
    GC-content and length.

-   Possibility to simulate PCR and sampling pseudo-replicates.

-   Simulation of size selection and polyadenylation with flexible
    target distributions.

-   Estimation of GC-dependent amplification efficiencies from real
    data, relying on assumptions about locality of biases and the mean
    efficiency of the fragment pool.

-   Estimation of relative expression levels.

-   Estimation of empirical fragment size distribution, model selection
    between normal vs. skew normal distributions.

-   Able to simulate experiments on the human transcriptome over a wide
    range of expression levels on a desktop machine.

Basic examples
--------------

The following basic examples can be run from the `src/` directory in the
package source tree:

-   Re-sample expression levels from a mixture of gamma distributions
    with two components with mean 5,000 and 10,000:
```
        $ ../tools/sel -d "0.5:g:(5000, 0.1) + 0.5:g:(10000, 100)" \
        test/basic/test_transcripts.fas > my_transcripts.fas
```

    The output my\_transcripts.fas is a Fasta file annotated with the
    expression levels:
```
        >ENST00000371588$9430
        GCTTCCGGCATCTGGCTCAGTTCCGCCATGGCCTCCTTGGAAGTCAGTCGTAGTCCTCGCAGGTCTCGGCGGGAGCTG ...
        ...
```

-   Simulate 10,000 fragments with default parameters, plot rlsim
    report:
```
        $ ./rlsim -n 10000 my_transcripts.fas > frags.fas
        $ ../tools/plot_rlsim_report
```

    The output file frags.fas contains the simulated fragments:
```
        >Frag_0 ENST00000374005 (Strand - Offset 1883 -- 2298)
        AGAGAATAGAGGGTAGAAGGGAAATTCTTGGCACCTGGACTAGAGTGAGATAAAAGGAGAGTAGGAAAGCAGTGA ...
        ...
```

-   Simulate paired-end sequencing using [simNGS](https://github.com/timmassingham/simNGS) and the
    runfile shipped with the package source:
```
        $ cat frags.fas | simNGS -p paired -o fastq -O reads test/cov/s_4_0066.runfile
```

    The files `reads_end1.fq` and `reads_end1.fq` contain the simulated
    paired-end reads:
```
        @Frag_24929 refB (Strand - Offset 11583 -- 12067) 151M
        GCCCCGAGTAGTTCTGGGCGGGGCCCCCGCGGCCAGCGCCGCCCACTATATATTATTTATTCTAACTATT ...
        +
        GGGEFG(EGFGGDGBFEGGG;FD7GEGGDGGA=GGGDG7FFGGGGGFGCAGGGGGFF?GGGFG@GGAGGG ...
        ...
```

-   Simulate 500,000 fragments, skew normal fragment size distribution
    with a spike, `after_prim_double` fragmentation method, 15 PCR
    cycles with the specified efficiency parameters, using 4 cores,
    verbose mode:
```
        $ ./rlsim -n 50000 -d "0.9:sn:(600,50,4,300,2000) + 0.1:n:(700,1,600,2000)" \
        -f after_prim_double -c 15 -eg "(1.0,0.5,0.8)" -el "(0.1,0.7,1.0)" \
        -t 4 -v my_transcripts.fas > frags.fas
        $ ../tools/plot_rlsim_report
```

-   Estimate parameters from SAM file sorted by read name using (verbose
    mode):
``
        $ ../tools/effest -v -f ../tools/test/ref.fas ../tools/test/aln1.sam 
```

-   Estimate parameters from SAM file sorted by read name using – assume
    15 PCR cycles and a pool efficiency of 0.9:
```
        $ ../tools/effest -v -c 15 -m 0.9 -f ../tools/test/ref.fas ../tools/test/aln1.sam 
```
-   **Simulate 20,000 fragments using the raw parameters estimated by `effest`,
    set minimum GC-dependent efficiency to 0.5 (verbose mode)**:
```
        $ ./rlsim -v -n 20000 -j raw_params.json -jm 0.5 my_transcripts.fas > frags.fas
        $ ../tools/plot_rlsim_report
```

Getting more help
-----------------

Please consult the [package documentation](http://bit.ly/rlsim-doc) for more help on the tools and the technical background.

The BioStar Q&A forum ([http://www.biostars.org](http://www.biostars.org)) is an excellent place to get additional help. 
The author of the package will monitor the posts having the `rlsim` tag.

Dependencies
------------

The package runs on 64-bit GNU/Linux operating systems. The rlsim tool
is written in [golang](http://golang.org) and shipped as a statically linked
executable for the `amd64` Linux platforms, hence it has no
dependencies.

The rlsim tool can be built for other architectures supported by the
compiler, however only the `amd64` architecture is supported and the
32-bit binaries might not work properly.

The parameter estimation tool (`effest`) and the additional tools are
written in `Python 2.x` and depend on a couple of packages:

-   [numpy](https://pypi.python.org/pypi/numpy) >= 1.6.2

-   [matplotlib](https://pypi.python.org/pypi/matplotlib) >= 1.1.0

-   [scipy](https://pypi.python.org/pypi/scipy) >= 0.10.1

-   [biopython](https://pypi.python.org/pypi/biopython) >= 1.60

-   A modified version of the [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html) package available
    from the GitHub repository under [https://github.com/sbotond/rlsim/tree/master/misc](https://github.com/sbotond/rlsim/tree/master/misc).

    **The official releases of the HTSeq package contain a bug causing
    segmenation fault when parsing certain paired-end datasets. Please
    use the modified version from the rlsim repository!**

These packages (with the exception of the modified HTSeq) are readily
installable from the [Python Package Index](https://pypi.python.org/pypi) using the [pip](https://pypi.python.org/pypi/pip)
tool by issuing the following command:

```
    pip install numpy matplotlib scipy biopython
```

The additional tool tools depend on samtools [samtools](http://samtools.sourceforge.net/). Simulating
Illumina sequencing of the fragments can be done by using the [simNGS](https://github.com/timmassingham/simNGS)
package (recommended) or any other sequencing simulator.


Installing the latest release
-----------------------------

The release tarballs can be obtained from the [releases](https://github.com/sbotond/rlsim/tree/master/releases) directory:

```
    wget https://github.com/sbotond/rlsim/blob/master/releases/rlsim-latest_amd64.tar.gz?raw=true t.tgz -O rlsim-latest_amd64.tar.gz
```

The unpacked release directory contains the following files:

-   `bin/` – directory containing the executables:

    -   `rlsim` – the main tool simulating library construction

    -   `effest` – tool for estimating selected parameters from real datasets

    -   `sel` – tool for sampling expression levels

    -   `plot_rlsim_report`– tool for plotting the `rlsim` report

    -   `pb_plot` – tool for visualising sequence biases

    -   `cov_cmp` – tool for comparing coverage trends across datasets

    -   `plot_cov` – tool for plotting transcript read coverage colored by
        reference base

-   COPYING – GPL v3 licence

-   README.md – short instructions in markdown format

-   `rlsim_manual.pdf` – package manual

The executables under `bin/` can be installed by copying them to a
directory listed in the `$PATH` environmental variable.

Building from source
--------------------

Build dependencies:

-   The package is built using `make` in a standard `Linux` environment.

-   Building the rlsim tool requires the `Go` compiler which can be installed
    as described on the projects [website](http://golang.org/doc/install).

-   A standard `LaTeX` installation with `pdflatex` has to be present in
    order to produce the package documentation.

-   Regenerating the test datasets needs the bwa [bwa](http://bio-bwa.sourceforge.net) and 
    [samtools](http://samtools.sourceforge.net) commands to be installed.

The package source can be obtained by cloning the `GitHub` [repository](https://github.com/sbotond/rlsim)
and built by issuing `make` in the top level directory:

```
    git clone https://github.com/sbotond/rlsim.git
    cd rlsim
    make
```

A tarball can be built by issuing:

```
   make release
```

## Experimental builds with gccgo

`rlsim` can be compiled with the recent version (>=4.7.2) of the `gccgo` compiler:

- Install the gc `Go` compiler suite on the projects [website](http://golang.org/doc/install), as the build process uses the `go` tool.

- Install gccgo as described [here](http://golang.org/doc/install/gccgo), or through your [package manager](http://packages.ubuntu.com/search?keywords=gccgo).

- Issue `make gccbuild` under `src/'.

[Benchmarks indicate](http://bit.ly/160TBUn) that the `gccgo` build is faster on average, however the difference in runtime is not substantial.

__Please note that the experimental `gccgo` builds are not supported. Feel free to use them, but please do not submit bug reports if anything goes wrong.__

Quick reference
---------------

### `rlsim`

```
    Simulate RNA-seq library preparation with priming biases, PCR biases and size selection (version: 1.3).

    Usage:
            rlsim [arguments] [transcriptome fasta files (optional)]

    Optional arguments:
                    argument                    type    default  
            -n      requested fragments         int     
            -d      fragment size distribution  string  "1.0:sn:(189, 24, -1.09975, 76, 294)" 
            -f      fragmentation method        string  "after_prim_double"
            -b      strand bias                 float   0.5
            -c      PCR cycles                  int     11
            -p      priming bias parameter      float   5.0
            -k      primer length               int     6
            -a      poly(A) tail size dist.     string  [check source]
            -flg    fragment loss probability   float   0.0
            -m      expression level multiplier float   1.0
            -e      fixed PCR efficiency        float   0.0
            -eg     GC efficiency parameters 
                    as "(shape, min, max)":     raw from SRR521457
            -el     length efficiency parameters 
                    as "(shape, min, max)":
                        shape                   float   0.0
                        min                     float   1.0
                        max                     float   1.0
            -j      raw parameter file          string  
                    superseeds -d, -c, -eg
            -jm     minimum raw gc efficiency   float   0.0
            -r      report file                 string  "rlsim_report.json"
            -t      number of cores to use      int     4
            -g      keep fragments in memory    bool    false
            -si     initial random seed         int     from UTC time
            -sp     pcr random seed             int     auto
            -ss     sampling random seed        int     auto
            -gobdir fragment directory          string  "rlsim_gob_$PID"
            -v      toggle verbose mode         bool    false
            -h      print usage and exit        bool    false
            -V      print version and exit      bool    false
            -prof   write CPU profiling info    string  ""
            -gcfreq trigger garbage collection  int     100
                    after this many transcripts
            -randt  generate RNG test files     bool    false

    Examples:
            rlsim -n 2000000 transcripts.fa
            cat transcripts.fa | rlsim -n 2000000

    For more details consult the package manual at:
            https://github.com/sbotond/rlsim/tree/master/doc/rlsim_manual.pdf
```

### `effest`

```
    usage: effest [-h] [-f ref_fasta] [-i iso_list] [-c nr_cycles] [-m mean_eff]
                  [-M max_eff] [-d dist_fam] [-g out_fasta] [-j out_json]
                  [-e expr_mul] [-a] [-t] [-w step_size] [-s out_count_file]
                  [-k in_count_file] [-p out_prior_file] [-o in_prior_file]
                  [-q min_qual] [-r report_file] [-l log_file] [-v]
                  [input file]

    Estimate GC-dependent fragment amplification efficiencies and fragment size
    distribution from paired-end RNA-seq data mapped to transcriptome (version
    1.1).

    positional arguments:
      input file         Aligned *paired end* reads in SAM format sorted by
                         *name*.

    optional arguments:
      -h, --help         show this help message and exit
      -f ref_fasta       Reference fasta.
      -i iso_list        List of single isoform transcripts.
      -c nr_cycles       Number of PCR cycles (11).
      -m mean_eff        Assumed pool efficiency (0.87).
      -M max_eff         Assumed maximum efficiency (None).
      -d dist_fam        Distribution to model fragment size distribution
                         (n|sn|*auto*).
      -g out_fasta       Output fasta.
      -j out_json        File to store estimated raw parameters (raw_params.json).
      -e expr_mul        Expression level multiplier (10000.0).
      -a                 Do not use GC efficiency correction on expression levels
                         (False).
      -t                 Trim off old expression values (True).
      -w step_size       Sliding window size / step size ratio (5).
      -s out_count_file  Pickle counts to the specified file (effest_counts.pk).
      -k in_count_file   Load counts from specifies file.
      -p out_prior_file  Pickle fragment prior to the specified file
                         (effest_pr.pk).
      -o in_prior_file   Load fragment prior from the specified pickle file.
      -q min_qual        Minimum mapping quality (0).
      -r report_file     Report PDF (effest_report.pdf).
      -l log_file        Log file.
      -v                 Toggle verbose mode (False).
```

### `plot_rlsim_report`

```
    usage: plot_rlsim_report [-h] [input file]

    Plot rlsim report (version 1.0).

    positional arguments:
      input file  rlsim report file.

    optional arguments:
      -h, --help  show this help message and exit
```

### `sel`

```
   usage: sel [-h] [-t] [-d dist_param] [-b nr_bins] [-r report_fil]
               [input fasta file]

    Sample expression levels from a mixture of gamma distributions (version 1.0).

    positional arguments:
      input fasta file  Transcripts and expression levels in Fasta format.

    optional arguments:
      -h, --help        show this help message and exit
      -t                Trim sequence names.
      -d dist_param     Expression level distribution.
      -b nr_bins        Number of bins in histogram.
      -r report_fil     Report PDF file.
```

### `pb_plot`

```
    usage: pb_plot [-h] [-f ref_fasta] [-r report_file] [-w winsize] [-i tr_list]
                   [-q min_qual] [-p pickle_file] [-s]
                   [input file]

    Visualise sequence biases around fragment start/end (version 1.1).

    positional arguments:
      input file      Aligned *paired end* reads in SAM format.

    optional arguments:
      -h, --help      show this help message and exit
      -f ref_fasta    Reference sequences in fasta format.
      -r report_file  Name of PDF report file.
      -w winsize      Window size.
      -i tr_list      List of single isoform transcripts.
      -q min_qual     Minimum mapping quality.
      -p pickle_file  Results pickle file.
      -s              Assume single ended dataset.
```

### `cov_cmp`

```
    usage: cov_cmp [-h] -f ref_fasta [-g] [-t nr_top] [-c min_cov] [-i iso_list]
                   [-l min_length] [-x] [-y] [-r report_file] [-q min_qual]
                   [-p pickle_file] [-v] [-s]
                   input file input file

    Compare relative coverage trends between the *expressed* transcripts of two
    datasets (version 1.1).

    positional arguments:
      input file      Two sets of aligned *paired end* reads in SAM format.

    optional arguments:
      -h, --help      show this help message and exit
      -f ref_fasta    Reference sequences in fasta format.
      -g              Do not color by AT/GC.
      -t nr_top       Plot at least this many top matches (30).
      -c min_cov      Minimum number of fragments per transcript (20).
      -i iso_list     List of single isoform genes.
      -l min_length   Minimum transcript length.
      -x              Sort by correlation coefficients.
      -y              Plot pairwise cumulative coverage.
      -r report_file  Name of PDF report file.
      -q min_qual     Minimum mapping quality (0).
      -p pickle_file  Results pickle file.
      -v              Toggle verbose mode.
      -s              Assume single ended dataset.

```

### `plot_cov`

```
    usage: plot_cov [-h] -r ref_fasta -b bam [-o outfile]

    Plot read coverage colored by the reference base (AT - blue, GC - red). This
    tools requires samtools to be installed in path.

    optional arguments:
      -h, --help    show this help message and exit
      -r ref_fasta  Reference transcriptome.
      -b bam        Position sorted and indexed BAM file.
      -o outfile    Output PDF (plot_cov.pdf).
```
            

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/a7401ff83389dacf5bf79399ba749f32 "githalytics.com")](http://githalytics.com/sbotond/rlsim)

