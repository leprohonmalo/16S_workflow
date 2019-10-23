# Script TP métagénomique 1
dossier_reads_bruts=$1
dossier_sortie=$2
dossier_init=$(pwd)


mkdir fastqc/
mkdir raw_data/
mkdir trimmed/
mkdir fasta/
mkdir dereplication/
mkdir chimera_filtering/
mkdir clustering/
mkdir OTU_table/


cd $dossier_reads_bruts

list_fastq=$(ls *fastq.gz)

cd $dossier_init

for i in $list_fastq; do
    echo $i
    cp $dossier_reads_bruts$i raw_data/$i
done 

for i in $list_fastq; do
    fastqc -o fastqc/ --noextract raw_data/$i
    gzip -d raw_data/$i
done

cd raw_data/
list_R1=$(ls *R1.fastq)
cd ../

for i in $list_R1; do
    echo fastq/$i
    echo fastq/"${i/R1.fastq/R2.fastq}"
    java -jar soft/AlienTrimmer.jar -if fastq/$i -ir fastq/"${i/R1.fastq/R2.fastq}" -of trimmed/$i -or trimmed/"${i/R1.fastq/R2.fastq}" -q 20 -c databases/contaminants.fasta
done

for i in $list_R1; do
    mv trimmed/$i trimmed/"${i/R1.fastq/R1.trimmed.fastq}"
    mv trimmed/"${i/R1.fastq/R2.fastq}" trimmed/"${i/R1.fastq/R2.trimmed.fastq}"
done

cd trimmed/
list_R1=$(ls *R1.trimmed.fastq)
cd ../

for i in $list_R1; do 
    echo ${i/R1.trimmed/R2.trimmed}
    echo ${i/R1.trimmed.fastq/trimmed.fasta}
    echo ";sample="${i/_R1.trimmed.fastq/}";"
    fastqc -o fastqc/ --noextract trimmed/$i
    fastqc -o fastqc/ --noextract trimmed/${i/R1.trimmed/R2.trimmed}
    vsearch --fastq_mergepairs trimmed/$i --reverse trimmed/${i/R1.trimmed/R2.trimmed} --fastaout fasta/${i/R1.trimmed.fastq/trimmed.fasta} --label_suffix ";sample="${i/_R1.trimmed.fastq/}";"
done

list_fasta=$(ls fasta/*)
echo $list_fasta

cat $list_fasta | tr -d "[:blank:]" > fasta/amplicon.fasta

vsearch --derep_fulllength fasta/amplicon.fasta --output dereplication/dereplication.fasta --sizeout --minuniquesize 10

vsearch --uchime_denovo dereplication/dereplication.fasta --nonchimeras chimera_filtering/non_chimeras.fasta --chimeras chimera_filtering/chimeras.fasta

vsearch --cluster_size chimera_filtering/non_chimeras.fasta --centroids clustering/clustered_seq.fasta --id 0.97 --clusterout_id --consout clustering/consensus_seq.fasta --sizeout --relabel "OTU_"

vsearch --usearch_global fasta/amplicon.fasta --db clustering/clustered_seq.fasta --otutabout OTU_table/OTU_table.txt --id 0.97 --sizeout

vsearch --usearch_global clustering/clustered_seq.fasta --db databases/mock_16S_18S.fasta --id 0.9 --top_hits_only --userfields query+target --userout OTU_table/annot_table.txt

mv fastqc/ $2
mv raw_data/ $2
mv trimmed/ $2
mv fasta/ $2
mv dereplication/ $2
mv chimera_filtering/ $2
mv clustering/ $2
mv OTU_table/ $2
