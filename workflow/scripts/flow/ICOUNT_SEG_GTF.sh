#!/bin/bash -euo pipefail   
iCount-Mini segment \   
    Homo_sapiens.GRCh38.109_bracketsremoved.cmd.gtf \   
    Homo_sapiens_seg.gtf \   
    Homo_sapiens.GRCh38.fasta.fai   
   
mv regions.gtf.gz Homo_sapiens_regions.gtf.gz   
   
cat <<-END_VERSIONS > versions.yml   
"PREPARE_CLIPSEQ:ICOUNT_SEG_GTF":   
    iCount-Mini: $(iCount-Mini -v)   
END_VERSIONS   