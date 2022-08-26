
<img src="man/figures/simmethods_logo.png" align="right" width = "210px" height="180px"/>

# A collection of 42 simulation methods for single-cell RNA-seq data

Simmethods collects and documents 36 popular and common simulation
methods for single-cell transcriptomics data. To satisfy user’s
requirements in different scenarios, we bundled the simulators
comprehensively and users can simulate many kinds of single-cell RNA-seq
data (different number of groups, batches, differential expressed genes
and even data with differentiation trajectory) using certain methods. If
you want to learn how to use a certain simulation method, please check
the following chart and commit an issue or send an Email to us when you
have any question.

## The list of simulation methods

| method         | language                                                           | url                                                                                                                                                    | doi                                                                                                                              | journal                    |
|:---------------|:-------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------|:---------------------------|
| BASiCS         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/BASiCS.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>             | <a href='https://doi.org/10.1371/journal.pcbi.1004333'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>     | PLoS Computational Biology |
| BEARscc        | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://www.bioconductor.org/packages/release/bioc/html/BEARscc.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>        | <a href='https://doi.org/10.1038/s41467-018-03608-y'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Nature Communications      |
| CancerInSilico | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://www.bioconductor.org/packages/release/bioc/html/CancerInSilico.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a> | <a href='https://doi.org/10.1371/journal.pcbi.1006935'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>     | PLoS Computational Biology |
| dropsim        | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/marchinilab/dropsim'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                      |                                                                                                                                  | NA                         |
| dyngen         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://cran.r-project.org/web/packages/dyngen/index.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                   | <a href='https://doi.org/10.1038/s41467-021-24152-2'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Nature Communications      |
| dyntoy         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/dynverse/dyntoy'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                          |                                                                                                                                  | NA                         |
| ESCO           | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/JINJINT/ESCO'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                             | <a href='https://doi.org/10.1093/bioinformatics/btab116'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>   | Bioinformatics             |
| hierarchicell  | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/kdzimm/hierarchicell'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                     | <a href='https://doi.org/10.1186/s12864-021-07635-w'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | BMC Genomics               |
| Kersplat       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-017-1305-0'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| Lun            | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-017-1305-0'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| Lun2           | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-017-1305-0'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| MFA            | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/kieranrcampbell/mfa'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                      | <a href='https://doi.org/10.12688/wellcomeopenres.11087.1'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a> | Wellcome Open Research     |
| muscat         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/HelenaLC/muscat'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                          | <a href='https://doi.org/10.1038/s41467-020-19894-4'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Nature Communications      |
| phenopath      | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/phenopath.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>          | <a href='https://doi.org/10.1101/159913'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>                   | bioRxiv                    |
| POWSC          | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='http://www.bioconductor.org/packages/release/bioc/html/POWSC.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1093/bioinformatics/btaa607'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>   | Bioinformatics             |
| powsimR        | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/bvieth/powsimR'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                           | <a href='https://doi.org/10.1093/bioinformatics/btx435'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>    | Bioinformatics             |
| PROSSTT        | <img src='man/figures/python_logo.png' height='18px' width='54px'> | <a href='http://wwwuser.gwdg.de/~compbiol/prosstt/doc/'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                               | <a href='https://doi.org/10.1093/bioinformatics/btz078'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>    | Bioinformatics             |
| scDD           | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://www.bioconductor.org/packages/release/bioc/html/scDD.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-016-1077-y'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| scDesign       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/Vivianstats/scDesign'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                     | <a href='https://doi.org/10.1093/bioinformatics/btz321'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>    | Bioinformatics             |
| scDesign2      | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/JSB-UCLA/scDesign2'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                       | <a href='https://doi.org/10.1186/s13059-021-02367-2'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Genome Biology             |
| SCRIP          | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://cran.r-project.org/web/packages/SCRIP/index.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                    |                                                                                                                                  | NA                         |
| Simple         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-017-1305-0'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| SparseDC       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://cran.rstudio.com/web/packages/SparseDC/index.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                   | <a href='https://doi.org/10.1093/nar/gkx1113'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>              | Nucleic Acids Research     |
| SPARSim        | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://gitlab.com/sysbiobig/sparsim'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                        | <a href='https://doi.org/10.1093/bioinformatics/btz752'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>    | Bioinformatics             |
| Splat          | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-017-1305-0'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| SplatPop       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://bioconductor.org/packages/release/bioc/html/splatter.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>           | <a href='https://doi.org/10.1186/s13059-021-02546-1'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Genome Biology             |
| SPsimSeq       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://www.bioconductor.org/packages/release/bioc/html/SPsimSeq.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>       | <a href='https://doi.org/10.1093/bioinformatics/btaa105'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>   | Bioinformatics             |
| SymSim         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/YosefLab/SymSim'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                          | <a href='https://doi.org/10.1038/s41467-019-10500-w'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Nature Communications      |
| TedSim         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/Galaxeee/TedSim'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                          | <a href='https://doi.org/10.1093/nar/gkac235'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>              | Nucleic Acids Research     |
| VeloSim        | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/PeterZZQ/VeloSim'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                         | <a href='https://doi.org/10.1101/2021.01.11.426277'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | bioRxiv                    |
| zinbwave       | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='http://www.bioconductor.org/packages/release/bioc/html/zinbwave.html'><img src='man/figures/URL.png' height='18px' width = '18px'></a>        | <a href='https://doi.org/10.1038/s41467-017-02554-5'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>       | Nature Communications      |
| zinbwaveZinger | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/statOmics/zinbwaveZinger'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                 | <a href='https://doi.org/10.1186/s13059-018-1406-4'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |
| zingeR         | <img src='man/figures/R_logo.png' height='18px' width='18px'>      | <a href='https://github.com/statOmics/zingeR'><img src='man/figures/URL.png' height='18px' width = '18px'></a>                                         | <a href='https://doi.org/10.1186/s13059-018-1406-4'><img src='man/figures/doi_logo.png' height='18px' width = '18px'></a>        | Genome Biology             |

## New methods

We are glad to add new simulation methods if some methods are innovative
and creative that many users commonly used. If you have the
requirements, please tell me by email (<duohongrui@cqnu.edu.cn>) or
raise an issue for that.
