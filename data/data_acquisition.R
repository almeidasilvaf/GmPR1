
#----Download soybean proteome----
library(here)
download.file("ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/Fasta/proteome.selected_transcript.gma.fasta.gz", 
              destfile = here("data", "Gmax_proteome.fa.gz"))

#----Download .gff file----
# Download gff from PLAZA
download.file(url="ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_04/GFF/gma/annotation.selected_transcript.exon_features.gma.gff3.gz", 
              destfile = here("data", "PLAZA_selected.transcripts.gff.gz"))
