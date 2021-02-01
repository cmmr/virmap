#!/usr/bin/env perl


use warnings;
use strict;

#why not?
my $bell = chr(7);

#make sure we aren't root
if ($< == 0) {
	die "don't run as root, this will probably work if you run as root, but it's not tested like that.";
}

#get amount of processors
my $procs = `nproc`;
chomp $procs;

#what linux am I?
my $linux;
my $user;
my $linuxFlavor = `cat /etc/*release* | grep ID_LIKE=`;
chomp $linuxFlavor;
if ($linuxFlavor =~ m/debian/ or $linuxFlavor =~ m/ubuntu/) {
	$linux = "deb";
	$user = "ubuntu";
} elsif ($linuxFlavor =~ m/centos/ or $linuxFlavor =~ m/rhel/) {
	$linux = "rh";
	$user = "ec2-user";
}
my $linuxID = `cat /etc/*release* | grep "^ID=" | sed "s/ID=//g" | tr -d '"'`;
chomp $linuxID;

unless ($linuxID and $linux) {
	die "unknown linux\n";
}

#update base machine
if ($linux eq "rh") {
	system("sudo yum -y update");
	if ($linuxID eq "amzn") {
		system("sudo yum install -y amazon-linux-extras");
		system("sudo amazon-linux-extras enable epel");
		system("sudo yum clean metadata");
	}
	system("sudo yum install -y epel-release");
	system("sudo yum -y groupinstall \"Development Tools\"");
	system("sudo yum -y install awscli cpan openssl-devel snappy snappy-devel zlib zlib-devel bzip2 bzip2-devel lz4-devel zstd.x86_64 lbzip2 lftp parallel java-devel pigz python3 nvme-cli expat-devel");
	system("sudo yum -y update");
	system("sudo pip3 install --upgrade pip");
} elsif ($linux eq "deb") {
	system("sudo apt-get update");
	system("sudo apt-get -y upgrade");
	system("sudo apt-get -y install build-essential");
	system("sudo apt-get -y install zstd libzstd-dev liblz4-dev libsnappy-dev zlib1g-dev libbz2-dev lbzip2 lftp parallel openjdk-8-jdk-headless pigz python3 python3-pip awscli nvme-cli libssl-dev libexpat1-dev");
	system("sudo apt-get update");
	system("sudo apt-get -y upgrade");
}


#make raid-0 scratch disk
my $nvmeList = `sudo nvme list | grep "Amazon EC2 NVMe Instance Storage" | cut -f1 -d " " | tr "\n" " "`;
chomp $nvmeList;
my $nvmeCount = `sudo nvme list | grep "Amazon EC2 NVMe Instance Storage" | wc -l`;
chomp $nvmeCount;
system("sudo mkdir /scratch");
if ($nvmeCount == 0) {
	die "no attached NVMe instance storage\n";
} elsif ($nvmeCount == 1) {
	system("sudo mkfs.ext4 $nvmeList");
	system("sudo mount $nvmeList /scratch");
} else {	
	system("sudo mdadm --create /dev/md0 --level=0 --raid-devices=$nvmeCount $nvmeList");
	system("sudo mkfs.ext4 /dev/md0");
	system("sudo mount /dev/md0 /scratch");
}
system("sudo chown -R $user:$user /scratch");
system("sudo chmod -R 0775 /scratch");
system("mkdir /scratch/tmp");
system("sudo chown -R $user:$user /usr/local");
system("sudo chmod -R 0775 /usr/local");


#set environment variables
$ENV{'PERL_MM_USE_DEFAULT'} = 1;
$ENV{'PATH'} = "/usr/local/bin:" .  $ENV{'PATH'};
$ENV{'TMPDIR'} = "/scratch/tmp";
if ($linux eq "deb") {
	$ENV{'JAVA_HOME'} = "/usr/lib/jvm/java-8-openjdk-amd64";
	$ENV{'PATH'} = "$ENV{'JAVA_HOME'}/bin:$ENV{'PATH'}";
	open OUT, ">>/home/$user/.bash_profile";
	print OUT "export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64\n";
	print OUT "export PATH=\$JAVA_HOME/bin:\$PATH";
	close OUT;
}
my $tmpdir = $ENV{'TMPDIR'};

#update coreutils
system("wget ftp://ftp.gnu.org/gnu/coreutils/coreutils-8.32.tar.gz -O $tmpdir/coreutils-8.32.tar.gz");
system("tar zxvf $tmpdir/coreutils-8.32.tar.gz -C $tmpdir");
system("cd $tmpdir/coreutils-8.32 && ./configure --prefix=/usr/local");
system("cd $tmpdir/coreutils-8.32 && make -j$procs");
system("cd $tmpdir/coreutils-8.32 && make install");

#update perl 
system("wget https://www.cpan.org/src/5.0/perl-5.30.2.tar.gz -O $tmpdir/perl-5.30.2.tar.gz");
system("tar zxvf $tmpdir/perl-5.30.2.tar.gz -C $tmpdir");
system("cd $tmpdir/perl-5.30.2 && ./Configure -Dcc=gcc -Dusethreads -Duselargefiles -Dusemorebits -des");
system("cd $tmpdir/perl-5.30.2 && make -j$procs");
system("cd $tmpdir/perl-5.30.2 && make test");
system("cd $tmpdir/perl-5.30.2 && make install");

#install required perl packages
system("cpan CPAN");
system("cpan App::cpanminus");
system("cpanm App::cpanoutdated");
system("cpan-outdated -p | cpanm");
system("cpanm Compress::Zstd");
system("cpanm REST::Client");
system("cpanm OpenSourceOrg::API");
system("cpanm POSIX::1003::Sysconf");
system("cpanm Sereal");
system("cpanm Text::Levenshtein::Damerau::XS");
system("cpanm Text::Levenshtein::XS");
system("cpanm Array::Shuffle");
system("cpanm Array::Split");
system("cpanm Sys::MemInfo");
system("cpanm XML::Hash::XS");
system("cpanm Cpanel::JSON::XS");
system("cpanm Statistics::Basic");
system("cpanm XML::Simple");


#install other software
system("mkdir -p /usr/local/virmap/src/");
system("git clone https://github.com/cmmr/virmap.git /usr/local/virmap/src/virmap");
system("ls /usr/local/virmap/src/virmap/DatabaseConstructor/*.pl | xargs -I {} ln -s {} /usr/local/bin/");
system("ls /usr/local/virmap/src/virmap/*.pl | xargs -I {} ln -s {} /usr/local/bin/");
system("wget http://korflab.ucdavis.edu/Unix_and_Perl/FAlite.pm -O /usr/local/lib/perl5/site_perl/5.30.2/FAlite.pm");
system("wget https://github.com/bbuchfink/diamond/releases/download/v0.9.32/diamond-linux64.tar.gz -O $tmpdir/diamond-linux64.tar.gz");
system("tar zxvf $tmpdir/diamond-linux64.tar.gz -C /usr/local/bin/");
system("tar zxvf /usr/local/virmap/src/virmap/Installer/BBMap_38.84.tar.gz -C /usr/local/virmap/src/");
system("cd /usr/local/virmap/src/bbmap/jni; make -f makefile.linux");
system("sudo cp /usr/local/virmap/src/bbmap/jni/libbbtoolsjni.so  /usr/lib/");
system("ls /usr/local/virmap/src/bbmap/*.sh | xargs -I {} ln -s {} /usr/local/bin/");
system("pip install wheel");
system("pip install khmer");
system("wget https://github.com/torognes/vsearch/releases/download/v2.14.2/vsearch-2.14.2-linux-x86_64.tar.gz -O $tmpdir/vsearch-2.14.2-linux-x86_64.tar.gz");
system("tar zxvf $tmpdir/vsearch-2.14.2-linux-x86_64.tar.gz -C /usr/local/virmap/src/");
system("ln -s /usr/local/virmap/src/vsearch-2.14.2-linux-x86_64/bin/vsearch /usr/local/bin/");
system("wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz -O $tmpdir/ncbi-blast-2.10.0+-x64-linux.tar.gz");
system("tar zxvf $tmpdir/ncbi-blast-2.10.0+-x64-linux.tar.gz -C /usr/local/virmap/src/");
system("ls /usr/local/virmap/src/ncbi-blast-2.10.0+/bin/* | xargs -I {} ln -s {} /usr/local/bin/");
system("wget https://github.com/voutcn/megahit/releases/download/v1.1.4/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz -O $tmpdir/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz");
system("tar zxvf $tmpdir/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz -C /usr/local/virmap/src/");
system("ls /usr/local/virmap/src/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin/megahit* | xargs -I {} ln -s {} /usr/local/bin/");
system("wget ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz -O $tmpdir/edirect.tar.gz");
system("tar zxvf $tmpdir/edirect.tar.gz -C /usr/local/virmap/src/");
system("ln -s /usr/local/virmap/src/edirect/efetch /usr/local/bin/");
system("ln -s /usr/local/virmap/src/edirect/edirect.pl /usr/local/bin/");




if ($linux eq "rh") {
	system("wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.6/setup-yum.sh -O $tmpdir/setup-yum.sh");
	system("chmod 0775 $tmpdir/setup-yum.sh");
	system("sudo $tmpdir/setup-yum.sh");
} elsif ($linux eq "deb") {
	system("wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.6/setup-apt.sh -O $tmpdir/setup-apt.sh");
	system("chmod 0775 $tmpdir/setup-apt.sh");
	system("sudo $tmpdir/setup-apt.sh");
}
$ENV{'PATH'} = "/usr/local/ncbi/sra-tools/bin:$ENV{'PATH'}";
open OUT, ">>/home/$user/.bash_profile";
print OUT "source /etc/profile.d/sra-tools.sh";
close OUT;


#install RocksDB for perl
system("wget http://www.cpan.org/authors/id/J/JI/JIRO/RocksDB-0.05.tar.gz -O $tmpdir/RocksDB-0.05.tar.gz");
system("tar zxvf $tmpdir/RocksDB-0.05.tar.gz -C $tmpdir");
system("patch $tmpdir/RocksDB-0.05/builder/MyBuilder.pm -i /usr/local/virmap/src/virmap/Installer/MyBuilder.pm.patch");
system("cd $tmpdir/RocksDB-0.05/; perl Build.PL; perl Build install");


bells(1);


sub bells {
	my $amount = $_[0];
	for (my $i = 0; $i < $amount; $i++) {
		print STDERR "$bell";
		sleep(1);
	}

}



__END__














