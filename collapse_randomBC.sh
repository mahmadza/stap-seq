#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

# Programs:
# * bowtie
# * slippage_filter_se.sh
# * Jim Kent's utils

# Files:
# * /groups/stark/indices/bowtie/${ASSEMBLY}/${ASSEMBLY}
# * /groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes

################################################################################
# Set default values
################################################################################

ASSEMBLY="dm3"
RAN_BARCODE_LEN=8          #random barcode length

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
    echo >&2 "
$(basename $0)  - collapse mapped fragments on positions and random barcodes
                - will also remove fragments with N in their barcodes

USAGE: $(basename $0) -i <BigBed input file> -o <BigBed output file> [OPTIONS]
 -i     Input file (bigBed) [required]
        Note: the mapped paired-end reads should be in column 4, format: FORWARDREAD_REVERSEREAD
        the random barcode will be extracted from 5' of the FORWARDREAD
 -o     Output file (bigBed) [required]
 -R     length of random barcode [default: $RAN_BARCODE_LEN]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
"
    exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "i:o:R:g:" o
do
    case "$o" in
        i) INFILE="$OPTARG";;
        o) OUTFILE="$OPTARG";;
        R) RAN_BARCODE_LEN="$OPTARG";;
        g) ASSEMBLY="$OPTARG";;
       \?) exit 1;;
    esac
done

if [ -z "$INFILE" -o -z "$OUTFILE" ]; then
    echo >&2 "ERROR: -i -o are required!"
    exit 1
fi

################################################################################
# Set index and chromosome sizes
################################################################################

if [ "$ASSEMBLY" = "dm3" ]; then
    SIZES=/groups/stark/genomes/chrom/dm3.chrom.sizes
else
    SIZES=/groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes
fi

[ -e "$SIZES" ] || echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"

frags_clean=$(mktemp)
bigBedToBed $INFILE stdout | \
    awk -vOFS="\t" '{
    random_bc=substr($4,1,8)
    if( random_bc ~ "N" ) next
    print $1,$2,$3,random_bc,$6
    }' | sort -k1,1 -k2,2n -k3,3n -k5,5 -k4,4 | \
        uniq -c | sort -k2,2 -k3,3n -k4,4n -k6,6 -k1,1nr | \
        awk -vOFS="\t" \
        'function bc_diff(bc_one,bc_two,schism){
            schism=0
            for(i=1;i<=length(bc_one);i++)
            {
                if(substr(bc_one,i,1)!=substr(bc_two,i,1))
                    schism++
            }
            return schism
        }
        
        (NR==1){
            last_position=$2" "$3" "$4" "$6
            old_barcode[$5]=$5
            print
            next
            }
            
        {
            if($2" "$3" "$4" "$6!=last_position)
            {
                #if encounters new position
    
                delete old_barcode
                last_position=$2" "$3" "$4" "$6
                old_barcode[$5]=1
                print
                next
            }
    
            else
    
            {
                #if not a new position...
    
                #go over all of the valid barcodes
                for(bc in old_barcode)
                    if(bc_diff($5,bc)==1)                   #skip if difference of two bc just 1
                        next
    
                #keep the valid barcode if difference >1
                old_barcode[$5]=1
                print
                
            }
        }' | awk -vOFS="\t" '{print $2,$3,$4,$5,0,$6,$1}' > $frags_clean

#convert to bigbed
bedToBigBed $frags_clean /groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes $OUTFILE

################################################################################
# Exit
################################################################################

rm -rf $frags_clean

exit 0

