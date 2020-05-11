# VirMAP

The full VirMAP source code is being cleaned and documented. It will be released as soon as possible. 

The paper is published at Nature Communications:
https://doi.org/10.1038/s41467-018-05658-8<br />
<br />

<b>Maximal viral information recovery from sequence data using VirMAP</b>

Nadim J Ajami1,2, Ω, Matthew C. Wong1,2, Ω, Matthew C. Ross1,2, Richard E. Lloyd2, Joseph F. Petrosino1,2.
 
1 Alkek Center for Metagenomics and Microbiome Research, and 2­ Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, Texas.
Ω ­Authors contributed equally to this work.
 
Corresponding author:
 
Nadim J. Ajami<br />
nadimajami@gmail.com
 
Abstract
 
Accurate classification of the human virome is critical to a full understanding of the role viruses play in health and disease. This implies the need for sensitive, specific, and practical pipelines that return precise outputs while still enabling case-specific post hoc analysis. Viral taxonomic characterization from metagenomic data suffers from high background noise and signal crosstalk that confounds current methods. Here we develop VirMAP that overcomes these limitations using techniques that merge nucleotide and protein information to taxonomically classify viral reconstructions independent of genome coverage or read overlap. We validate VirMAP using published data sets and viral mock communities containing RNA and DNA viruses and bacteriophages. VirMAP offers opportunities to enhance metagenomic studies seeking to define virome-host interactions, improve biosurveillance capabilities, and strengthen molecular epidemiology reporting.

Cite as:
Ajami, N. J., Wong, M. C., Ross, M. C., Lloyd, R. E., & Petrosino, J. F. (2018). Maximal viral information recovery from sequence data using VirMAP. Nature Communications, 9(1), 3205. https://doi.org/10.1038/s41467-018-05658-8
<br />

Rough Requirements:<br />

Hardware:<br />
48GB+ RAM<br />
4 cores+ CPU<br />

Perl:<br />
Multi-threaded 5.24+<br />

Installed programs:
diamond (https://github.com/bbuchfink/diamond)<br />
bbtools (https://jgi.doe.gov/data-and-tools/bbtools)<br />
blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST)<br />
lbzip2 (http://lbzip2.org/download)<br />
MEGAHIT (https://github.com/voutcn/megahit)<br />
normalize-by-median.py (https://github.com/dib-lab/khmer)<br />
pigz (https://github.com/madler/pigz)<br />
vsearch (https://github.com/torognes/vsearch)<br />
zstd (https://github.com/facebook/zstd)<br />

CPAN Dependencies:<br />
Compress::Zstd (https://metacpan.org/pod/Compress::Zstd)<br />
OpenSourceOrg::API (https://metacpan.org/pod/OpenSourceOrg::API)<br />
POSIX::1003::Sysconf (https://metacpan.org/pod/distribution/POSIX-1003/lib/POSIX/1003/Sysconf.pod)<br />
POSIX::RT::Semaphore (https://metacpan.org/pod/POSIX::RT::Semaphore)<br />
RocksDB (https://metacpan.org/pod/RocksDB)<br />
Sereal (https://metacpan.org/pod/Sereal)<br />
Text::Levenshtein::Damerau::XS (https://metacpan.org/pod/Text::Levenshtein::Damerau::XS)<br />
Text::Levenshtein::XS (https://metacpan.org/pod/Text::Levenshtein::XS)<br />
 

Custom Dependencies:<br />
FAlite (http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm)
