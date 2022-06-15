#!/bin/bash

#EGA_ACCESS="FILE_WITH_EGA_ACCESS_CREDENTIALS"
EGA_ACCESS="/home/smith10/ega_access"
EGA_DATA="EGAD00010000644"

if ! [ -d "data" ]; then
	mkdir data
fi

if ! [ -d "data/CCLE/downloaded" ]; then
	mkdir -p data/CCLE/downloaded
fi
wget -P data/CCLE/downloaded -c https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_SNP.Arrays_2012-10-30_1.tar.gz
tar -xzf data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_1.tar.gz -C data/CCLE/downloaded/
rm data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_1.tar.gz

wget -P data/CCLE/downloaded -c https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_SNP.Arrays_2012-10-30_2.tar.gz
tar -xzf data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_2.tar.gz -C data/CCLE/downloaded/
rm data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_2.tar.gz

wget -P data/CCLE/downloaded -c https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_SNP.Arrays_2012-10-30_3.tar.gz
tar -xzf data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_3.tar.gz -C data/CCLE/downloaded/
rm data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_3.tar.gz

wget -P data/CCLE/downloaded -c https://data.broadinstitute.org/ccle_legacy_data/dna_copy_number/CCLE_SNP.Arrays_2012-10-30_4.tar.gz
tar -xzf data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_4.tar.gz -C data/CCLE/downloaded/
rm data/CCLE/downloaded/CCLE_SNP.Arrays_2012-10-30_4.tar.gz

if ! [ -d "data/CCLE/celfiles" ]; then
        mkdir -p data/CCLE/celfiles
fi

find data/CCLE/downloaded/ -type f -iname "*.cel" | xargs -I {} mv {} data/CCLE/celfiles/

if ! [ -d "data/GDSC/downloaded" ]; then
        mkdir -p data/GDSC/downloaded
fi

pyega3 -cf ${EGA_ACCESS} fetch ${EGA_DATA} --saveto data/GDSC/downloaded/

if ! [ -d "data/GDSC/celfiles" ]; then
        mkdir -p data/GDSC/celfiles
fi

find data/GDSC/downloaded/ -type f -iname "*.cel" | xargs -I {} mv {} data/GDSC/celfiles/

echo -e "cel_files" > data/cel_files
find ./ -type f -iname *.cel >> data/cel_files

if ! [ -d "data/ascat_results" ]; then
        mkdir -p data/ascat_results
fi
