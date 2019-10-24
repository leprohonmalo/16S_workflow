# Script TP métagénomique 1
dossier_reads_bruts=$1
dossier_sortie=$2
dossier_init=$(pwd)

mkdir $2
mkdir $2fastqc/
mkdir $2raw_data/
mkdir $2trimmed/
mkdir $2fasta/
mkdir $2dereplication/
mkdir $2chimera_filtering/
mkdir $2clustering/
mkdir $2OTU_table/


cd $dossier_reads_bruts

list_fastq=$(ls *fastq.gz)

cd $dossier_init

for i in $list_fastq; do
    cp $dossier_reads_bruts$i $2raw_data/$i
done 

for i in $list_fastq; do
    fastqc -o $2fastqc/ --noextract $2raw_data/$i
    gzip -d $2raw_data/$i
done

cd $2raw_data 
list_R1=$(ls *R1.fastq)
cd $dossier_init

cd soft/
./JarMaker.sh AlienTrimmer.java
cd $dossier_init

for i in $list_R1; do
    java -jar soft/AlienTrimmer.jar -if $2raw_data/$i -ir $2raw_data/"${i/R1.fastq/R2.fastq}" -of $2trimmed/"${i/R1.fastq/R1.trimmed.fastq}" -or $2trimmed/"${i/R1.fastq/R2.trimmed.fastq}" -q 20 -c databases/contaminants.fasta
done

cd $2trimmed/
list_R1=$(ls *R1.trimmed.fastq)
cd $dossier_init

for i in $list_R1; do 
    fastqc -o $2fastqc/ --noextract $2trimmed/$i
    fastqc -o $2fastqc/ --noextract $2trimmed/${i/R1.trimmed/R2.trimmed}
    vsearch --fastq_mergepairs $2trimmed/$i --reverse $2trimmed/${i/R1.trimmed/R2.trimmed} --fastaout $2fasta/${i/R1.trimmed.fastq/trimmed.fasta} --label_suffix ";sample="${i/_R1.trimmed.fastq/}";"
done

list_fasta=$(ls $2fasta/*)

cat $list_fasta | tr -d "[:blank:]" > $2fasta/amplicon.fasta

vsearch --derep_fulllength $2fasta/amplicon.fasta --output $2dereplication/dereplication.fasta --sizeout --minuniquesize 10

vsearch --uchime_denovo $2dereplication/dereplication.fasta --nonchimeras $2chimera_filtering/non_chimeras.fasta --chimeras $2chimera_filtering/chimeras.fasta

vsearch --cluster_size $2chimera_filtering/non_chimeras.fasta --centroids $2clustering/clustered_seq.fasta --id 0.97 --clusterout_id --consout $2clustering/consensus_seq.fasta --sizeout --relabel "OTU_"

vsearch --usearch_global $2fasta/amplicon.fasta --db $2clustering/clustered_seq.fasta --otutabout $2OTU_table/OTU_table.txt --id 0.97 --sizeout

vsearch --usearch_global $2clustering/clustered_seq.fasta --db databases/mock_16S_18S.fasta --id 0.9 --top_hits_only --userfields query+target --userout $2OTU_table/annot_table.txt
