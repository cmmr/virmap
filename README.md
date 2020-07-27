# VirMAP

Please post any questions, comments, or concerns. 

The paper is published at Nature Communications:
https://doi.org/10.1038/s41467-018-05658-8<br />
<br />

<b>Maximal viral information recovery from sequence data using VirMAP</b>

Nadim J Ajami1,2, Ω, Matthew C. Wong1,2, Ω, Matthew C. Ross1,2, Richard E. Lloyd2, Joseph F. Petrosino1,2.
 
1 Alkek Center for Metagenomics and Microbiome Research, and 2­ Department of Molecular Virology and Microbiology, Baylor College of Medicine, Houston, Texas.<br />
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

Hardware (recommanded):<br />
64GB+ RAM<br />
12 cores+ CPU<br />
per instance.

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
REST::Client (https://metacpan.org/pod/REST::Client)<br />
OpenSourceOrg::API (https://metacpan.org/pod/OpenSourceOrg::API)<br />
POSIX::1003::Sysconf (https://metacpan.org/pod/distribution/POSIX-1003/lib/POSIX/1003/Sysconf.pod)<br />
RocksDB (https://metacpan.org/pod/RocksDB)<br />
Sereal (https://metacpan.org/pod/Sereal)<br />
Text::Levenshtein::Damerau::XS (https://metacpan.org/pod/Text::Levenshtein::Damerau::XS)<br />
Text::Levenshtein::XS (https://metacpan.org/pod/Text::Levenshtein::XS)<br />
Array::Shuffle (https://metacpan.org/pod/Array::Shuffle <br />
Array::Split (https://metacpan.org/pod/Array::Split)<br />
Sys::MemInfo (https://metacpan.org/pod/Sys::MemInfo)<br />
XML::Hash::XS (https://metacpan.org/pod/XML::Hash::XS)<br />
Cpanel::JSON::XS (https://metacpan.org/pod/Cpanel::JSON::XS)<br />
Statistics::Basic (https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod)<br />


Custom Dependencies:<br />
FAlite (http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm)<br />


## Installation instructions for AWS

### Suggested Requirements
r5d.24xlarge (recommended) or c5d.24xlarge.<br />
Amazon Linux 2 (recommended) or Ubuntu.<br />
Root volume size >30Gb<br \>
Local instance store SSD with >500 GB memory and 64GB RAM. As Genbank expands, the minimum SSD and RAM requirements will expand as well.<br />
Instructions have only been tested on a fresh Amazon Linux 2 image.<br />

### Initial Machine Setup

Example command line input for Amazon Linux 2:<br />

#### Update machine
`sudo yum -y update`

#### Download and run installer
`wget https://raw.githubusercontent.com/cmmr/virmap/master/virmapInstall.pl`<br />
`chmod 0775 virmapInstall.pl`<br />

#### Install Virmap with dependancies
`./virmapInstall.pl`

#### Create DB (can take > 1h)
`mkdir /scratch/VirmapDb`<br />
`makeVirmapDb.pl --outputDir /scratch/VirmapDb`

#### Set up sra-tools configuration
`vdb-config -i`<br />
Exit out of vdb-config by hitting 'x'

#### Grab Viral mock community
`mkdir /home/$USER/VirmapTest`<br />
`fasterq-dump -t /dev/shm -e 4 -O /home/$USER/VirmapTest SRR9875293`

#### Set TMPDIR to somewhere on /scratch
`mkdir /scratch/tmp`<br />
`export TMPDIR=/scratch/tmp`

#### Test run viral mock community
`Virmap.pl --threads $(nproc) --readF /home/$USER/VirmapTest/SRR9875293_1.fastq --readR /home/$USER/VirmapTest/SRR9875293_2.fastq --useMegahit --useBbnorm --sampleName SRR9875293 --outputDir /home/$USER/VirmapTest/VirmapRun --taxaJson /scratch/VirmapDb/Taxonomy.virmap --virDmnd /scratch/VirmapDb/virDmnd.dmnd --virBbmap /scratch/VirmapDb/virBbmap --gbBlastn /scratch/VirmapDb/gbBlastn --gbBlastx /scratch/VirmapDb/gbBlastx.dmnd 2>/home/$USER/VirmapTest/VirmapRun.err`


### Saving the installation for future use

Save VirmapDB in an s3 location.<br />
Make a snapshot/image of the VM.<br />
When running Virmap, copy the database from S3 to the local instance store SSD for use. Using the database over gp2 is not recommended.<br />




### Instructions for use after initial set up
Launch image.<br />
Machine must have at least 64GB RAM per simultaneous instance of Virmap.<br />
Local instance SSD >500GB.<br />

1. ssh into machine.

2. Recreate local scratch space.

`sudo mkdir /scratch; nvmeList = $(sudo nvme list | grep "Amazon EC2 NVMe Instance Storage" | cut -f1 -d " " | tr "\n" " "); nvmeCount = $(sudo nvme list | grep "Amazon EC2 NVMe Instance Storage" | wc -l); if [ $nvmeCount -gt 1 ]; then sudo mdadm --create /dev/md0 --level=0 --raid-devices=$nvmeCount $nvmeList; sudo mkfs.ext4 /dev/md0; sudo mount /dev/md0 /scratch; else sudo mkfs.ext4 $nvmeList; sudo mount $nvmeList /scratch; fi`

3. Set permissions on /scratch.

`sudo chown -R <user>:<group> /scratch`<br />
`sudo chmod -R 0775 /scratch`

4. Copy VirmapDB to scratch from your s3 space.

5. Set TMPDIR to somewhere on /scratch.



