#!/bin/bash -euo pipefail   
iCount-Mini segment \   
    Homo_sapiens_filtered.gtf \   
    Homo_sapiens_filtered_seg.gtf \   
    Homo_sapiens.GRCh38.fasta.fai   
   
mv regions.gtf.gz Homo_sapiens_filtered_regions.gtf.gz   
   
cat <<-END_VERSIONS > versions.yml   
"PREPARE_CLIPSEQ:ICOUNT_SEG_FILTGTF":   
    iCount-Mini: $(iCount-Mini -v)   
END_VERSIONS   