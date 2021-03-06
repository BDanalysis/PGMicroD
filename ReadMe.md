-----------------------------------------------------------------------------------------------
Installation

##Basic Environment

Linux operation system with python3.x

##Third Party Tools 

[Bwa]

a. Download: 

wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download ./

b. Unzip the file:

$ tar -xvf bwa-0.7.17.tar.bz2

c. Add Bwa into profile

$ vim .bashrc

$ export Dir_Bwa/bwa-0.7.17:$PATH      #Dir_Bwa is the abosulte directory of Bwa 

$ source .bashrc

[SAMtools]

a. Download:

wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2/download ./

b. Unzip the file:

$ tar -xvf samtools-1.7.tar.bz2

c. Add SAMTools into profile:

$ vim .bashrc

$ export Dir_SAMtools/samtools-1.7:$PATH     #Dir_SAMtools is the abosulte directory of SAMTools

$ source .bashrc

[EMBOSS]

a. Download:

ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
   
b. Unzip and install:

$ tar -xvf EMBOSS-6.6.0.tar.gz

$ cd EMBOSS-6.6.0

$ ./configure

$ make

(Attention: if error occurs : “…, without plugin: x11”, execute the following command:

$ ./configure --without-x

$ make

) 

$ sudo make install
   
c. Check whether installation is successful:

$ wossname -auto |more

(Attention:

If command usage appears, then EMBOSS is installed successfully; 

If error occurs: “…, wossname: error while loading shared libraries: libnucleus.so.6: cannot open shared object 

file: No such file or directory”, you should input the following command:

$ sudo /sbin/ldconfig    

)
    
d. Add EMBOSS into profile:

$ vim .bashrc

$ export Dir_EMBOSS/EMBOSS-6.6.0:$PATH              #Dir_EMBOSS is the abosulte directory of EMBOSS

$ source .bashrc


##Extra Python Library

$ pip3 install pysam

$ pip3 install sklearn

$ pip3 install scipy



-----------------------------------------------------------------------------------------------
Usage of software

if the sequencing reads is SE, please use "Microbe-SE";
if the sequencing reads is PE, please use "Microbe-PE".

1) Microbe-SE usage:

$ cd Microbe-SE/

$ sh MicrobeRun.sh   [FastaDir]   [ReadFileDir]   [ResultFileDir]

#ReadFileDir is the absolute director of SE reads file, ResultFileDir is the absolute director of result

2) Microbe-PE usage:

$ cd Microbe-PE/

$ sh MicrobeRun.sh   [FastaDir]   [ReadFile1Dir]   [ReadFile2Dir]   [ResultFileDir]

#ReadFile1Dir and ReadFile2Dir is the absolute director of PE1 and PE2 reads file in PE respecitvely


-----------------------------------------------------------------------------------------------
Result Analysis

The result is composed of tow columns, the first column is microbe components, the second column is microbe abundance.


-----------------------------------------------------------------------------------------------
Demos Analysis

The PE sequencing sample can be found in 'software/Microbe-PE/example';

the SE sequencing sample can be found in 'software/Microbe-SE/example'.


