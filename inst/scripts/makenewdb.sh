#!/usr/bin/env bash
eval "$(conda shell.bash hook)"

# Do you have a SRST2 env? Make it, add biopython, download DB scripts, add python shebang to scripts, copy to bin
#sbatch -p NMLResearch -c 8 --mem=12G --wrap="conda create -y -n srst2 srst2 biopython"
#mkdir $CONDA_PREFIX/software
#cd $CONDA_PREFIX/software
#wget https://raw.githubusercontent.com/katholt/srst2/master/database_clustering/cdhit_to_csv.py
#wget https://raw.githubusercontent.com/katholt/srst2/master/database_clustering/csv_to_gene_db.py
#sed -i '1 i #!/usr/bin/env python' cdhit_to_csv.py
#sed -i '1 i #!/usr/bin/env python' csv_to_gene_db.py
#ln -s $CONDA_PREFIX/software/* $CONDA_PREFIX/bin/

echo "Making and moving into directory `date +'%Y%m%d'`"
mkdir `date +'%Y%m%d'`
cd `date +'%Y%m%d'`

echo "Activate the SRST2 env"
conda activate srst2

if [ $# -eq 0 ]
  then
    echo "No fasta supplied, getting the latest refseq plasmids"
    exit 1
    wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz
    gunzip plasmid.1.1.genomic.fna.gz    
  else
     echo "Hello world"
     exit 1
fi

#Do you have a CDHIT environ?
# sbatch -p NMLResearch -c 1 --mem 1G --wrap="conda create -y -n CDHIT-maxseq10mn cd-hit=4.8.1"
#Yes? Great. Now install CDHIT from source within it. 
#mkdir $CONDA_PREFIX/software
#cd $CONDA_PREFIX/software/
#srun wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
#tar -xvf cd-hit-v4.8.1-2019-0228.tar --gunzip
#cd cd-hit-v4.8.1-2019-0228
#srun make MAX_SEQ=10000000
#cd cd-hit-auxtools/
#make MAX_SEQ=10000000
#mv $CONDA_PREFIX/bin/cd-hit-est $CONDA_PREFIX/bin/BKcd-hit-est
#ln -s $CONDA_PREFIX/software/cd-hit-v4.8.1-2019-0228/cd-hit-est $CONDA_PREFIX/bin/cd-hit-est

# Go back to your new db dir and cluster: Refseq
echo "Cluster at 99% length and identity using custom built cd-hit with MAXSEQ=10mn"
conda activate CDHIT-maxseq10mn
sbatch -p NMLResearch -c 32 --mem=64G --wrap="cd-hit-est -i plasmid.1.1.genomic.fna -o db.cluster -c 0.99 -n 10 -s 0.99 -M 0 -T 0"

echo "Waiting on CDHIT"
while [ ! -f db.cluster.clstr ]; do sleep 5; done

echo "Collapse the clusters into their representatives"
cp db.cluster.clstr db.collapsed.clstr
sed -i '/^[^>].*[^\*\s]$/d' db.collapsed.clstr

echo "Prepare SRST2 formatted databases using scripts"
conda activate srst2
mkdir fasta
cd fasta
sbatch -p NMLResearch -c 4 --mem=16G --wrap="cdhit_to_csv.py --cluster_file ../db.collapsed.clstr --infasta_file ../plasmid.1.1.genomic.fna --outfile ../db.csv"
while [ ! -f ../db.csv ]; do sleep 5; done

cat *.fsa > collapsed.fasta
sbatch -p NMLResearch -c 4 --mem=16G --wrap="csv_to_gene_db.py -t ../db.csv -o ../SRST2_`date +'%Y%m%d'`.fasta -c 3 -f collapsed.fasta"

echo "Waiting on Gene DB"
while [ ! -f ../SRST2_`date +'%Y%m%d'`.fasta ]; do sleep 5; done

echo "New Database SRST2_`date +'%Y%m%d'`.fasta is ready"