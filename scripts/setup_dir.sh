#!/bin/bash
APT="https://downloads.thermofisher.com/APT/APT_2.11.4/apt_2.11.4_linux_64_bit_x86_binaries.zip"
ASCAT="https://github.com/VanLoo-lab/ascat.git"
APT_LIBS="http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip"
PENNCNV="https://github.com/WGLab/PennCNV.git"
ASCATSC="https://github.com/VanLoo-lab/ASCAT.sc.git"

APT_BASE=$(basename ${APT})
wget -c -P "resources/" ${APT}
unzip -d resources/ resources/${APT_BASE}
rm resources/${APT_BASE}
chmod +rx resources/apt_2.11.4_linux_64_bit_x86_binaries/bin/*

git clone ${ASCAT}
mv ascat resources/

git clone ${ASCATSC} 
mv ASCAT.sc resources/

if ! [ -d "refs" ]; then
	mkdir refs
fi
if ! [ -d "logs" ]; then
        mkdir logs
fi

APT_LIBS_BASE=$(basename ${APT_LIBS})
wget -c -P "refs/" ${APT_LIBS}
unzip -d refs/ refs/${APT_LIBS_BASE}
rm refs/${APT_LIBS_BASE}

git clone ${PENNCNV}
mv PennCNV refs/
cp refs/PennCNV/affy/libgw6/* refs/CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/
cp refs/PennCNV/affy/bin/normalize_affy_geno_cluster.pl resources/apt_2.11.4_linux_64_bit_x86_binaries/bin/
cp refs/PennCNV/kcolumn.pl resources/apt_2.11.4_linux_64_bit_x86_binaries/bin/

Rscript -e 'devtools::install("resources/ascat/ASCAT")'
Rscript -e 'devtools::install("resources/ASCAT.sc/")'

