# VirMAP

The full VirMAP source code is being cleaned and documented. It will be released as soon as possible. 

The paper is published at Nature Communications:
https://doi.org/10.1038/s41467-018-05658-8

Nadim J Ajami1,2, Ω, Matthew C. Wong1,2, Ω, Matthew C. Ross1,2, Richard E. Lloyd2, Joseph F. Petrosino1,2.
 
1 Alkek Center for Metagenomics and Microbiome Research, and 2­ Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, Texas.
Ω ­Authors contributed equally to this work.
 
Corresponding author:
 
Nadim J. Ajami
nadimajami@gmail.com
 
Abstract
 
Accurate classification of the human virome is critical to a full understanding of the role viruses play in health and disease. This implies the need for sensitive, specific, and practical pipelines that return precise outputs while still enabling case-specific post hoc analysis. Viral taxonomic characterization from metagenomic data suffers from high background noise and signal crosstalk that confounds current methods. Here we develop VirMAP that overcomes these limitations using techniques that merge nucleotide and protein information to taxonomically classify viral reconstructions independent of genome coverage or read overlap. We validate VirMAP using published data sets and viral mock communities containing RNA and DNA viruses and bacteriophages. VirMAP offers opportunities to enhance metagenomic studies seeking to define virome-host interactions, improve biosurveillance capabilities, and strengthen molecular epidemiology reporting.


Rough Requirements:
48GB+ RAM
4 cores+
Multi-threaded perl 5.24+


CPAN Dependencies:
Compress::Zstd
OpenSourceOrg::API
POSIX::1003::Sysconf
POSIX::RT::Semaphore
RocksDB
Sereal
Text::Levenshtein::Damerau::XS
Text::Levenshtein::XS
 

Custom Dependencies:
FAlite (http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm)
