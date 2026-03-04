#!/bin/bash

## Setup Conda environment
echo "Creating conda environment..."
conda create --solver libmamba -y -n CLEVER
conda activate CLEVER
conda install -y anaconda::wget
conda install -y bioconda::prodigal
conda install -y bioconda::cd-hit
conda install -y bioconda::diamond
conda install -y bioconda::blast
mkdir dbs

## To download the files to the right locations, change these variables for where to find the relevant files:
RESFINDER='https://bitbucket.org/genomicepidemiology/resfinder_db/raw/eecf0aa207594fe6d51badf808473de62b28cb06/'
CARD='https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2'
RESFINDERFG='https://raw.githubusercontent.com/RemiGSC/ResFinder_FG_Construction/refs/heads/main/output/RFG_db/ResFinder_FG_AA.faa'
PLSDB='https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta'
VERSION=1
CLEVER_ANNOTATION='clever_annotation_2026-03-01.txt'

## Prepare ResFinder files
echo "Preparing ResFinder..."
mkdir dbs/ResFinder
wget $RESFINDER\all.fsa
wget $RESFINDER\VERSION
wget $RESFINDER\phenotypes.txt
mv all.fsa dbs/ResFinder/
mv VERSION dbs/ResFinder/
mv phenotypes.txt dbs/ResFinder/
prodigal -q -a dbs/ResFinder/all.faa -i dbs/ResFinder/all.fsa >> progigal_log.txt
sed -i "s/\*//g" dbs/ResFinder/all.faa
sed -i "s/^>/>ResFinder-/" dbs/ResFinder/all.faa

## Prepare CARD files
echo "Preparing CARD..."
mkdir dbs/CARD
wget $CARD
mv *.tar.bz2 dbs/CARD/
bunzip2 dbs/CARD/*.tar.bz2
tar -xvf dbs/CARD/*.tar
mv aro_* dbs/CARD/
mv card.json dbs/CARD/
mv nucleotide_* dbs/CARD/
mv protein_* dbs/CARD/
mv shortname_* dbs/CARD/
mv snps.txt dbs/CARD/
mv build.txt dbs/CARD/
mv CARD-Download-README.txt dbs/CARD/
mv PMID.tsv dbs/CARD/
cat dbs/CARD/protein_fasta_protein_knockout_model.fasta dbs/CARD/protein_fasta_protein_overexpression_model.fasta dbs/CARD/protein_fasta_protein_variant_model.fasta > dbs/CARD/not_ARGs.fasta
sed -i "s/^>/>CARD-/" dbs/CARD/protein_fasta_protein_homolog_model.fasta

## Prepare ResFinderFG files
echo "Preparing ResFinderFG..."
mkdir dbs/ResFinderFG
wget $RESFINDERFG
mv ResFinder_FG_AA.faa dbs/ResFinderFG/
sed -i "s/^>/>ResFinderFG-/" dbs/ResFinderFG/ResFinder_FG_AA.faa

## Prepare fARGene predictions from Inda-Diaz et al 2023
echo "Preparing other sources..."
echo "   - Inda-Diaz 2023"
prodigal -q -a dbs/Other_sources/Inda-Diaz_2023.faa -i dbs/Other_sources/Inda-Diaz_2023.fasta >> progigal_log.txt
sed -i "s/\*//g" dbs/Other_sources/Inda-Diaz_2023.faa
sed -i "s/^>/>Inda-Diaz_2023-/" dbs/Other_sources/Inda-Diaz_2023.faa

## Merge all other sources into one file
echo "   merging..."
cat dbs/Other_sources/*.faa > Other_sources.faa

## Prepare PLSDB files
echo "Preparing PLSDB..."
mkdir dbs/PLSDB
wget $PLSDB
mv download_fasta dbs/PLSDB/plsdb.fna

## Match and cluster VARIANT database (100% identity)
echo "Creating variant database..."
diamond makedb --in dbs/CLEVER/CLEVER.variants.faa --db CLEVER --ignore-warnings
diamond blastp -q dbs/ResFinder/all.faa --db CLEVER -o ResFinder_vs_CLEVER.blastp --id 100 --query-cover 99 --subject-cover 99 --outfmt 6 --un ResFinder_unique.faa
cat dbs/CLEVER/CLEVER.variants.faa ResFinder_unique.faa > ResFinder.faa
diamond makedb --in ResFinder.faa --db ResFinder --ignore-warnings
diamond blastp -q dbs/CARD/protein_fasta_protein_homolog_model.fasta --db ResFinder -o CARD_vs_ResFinder.blastp --id 100 --query-cover 99 --subject-cover 99 --outfmt 6 --un CARD_unique.faa
diamond makedb --in dbs/CARD/not_ARGs.fasta --db CARD_NOT_ARGs --ignore-warnings
diamond blastp -q CARD_unique.faa --db CARD_NOT_ARGs -o CARD_vs_CARD_NON_ARGs.blastp --id 70 --query-cover 80 --subject-cover 80 --outfmt 6 --un CARD_filtered.faa
cat dbs/ResFinder/all.faa CARD_filtered.faa > ResFinder+CARD.faa
diamond makedb --in ResFinder+CARD.faa --db ResFinder+CARD --ignore-warnings
diamond blastp -q dbs/ResFinderFG/ResFinder_FG_AA.faa --db ResFinder+CARD -o ResFinderFG_vs_ResFinder+CARD.blastp  --id 100 --query-cover 99 --subject-cover 99 --outfmt 6 --un ResFinderFG_unique.faa
cat dbs/ResFinder/all.faa CARD_filtered.faa ResFinderFG_unique.faa > ResFinder+CARD+FG.faa
diamond makedb --in ResFinder+CARD+FG.faa --db ResFinder+CARD+FG --ignore-warnings
diamond blastp -q Other_sources.faa --db ResFinder+CARD+FG -o OS_vs_ResFinder+CARD+OS.blastp --id 100 --query-cover 99 --subject-cover 99 --outfmt 6 --un OS_unique.faa
cat dbs/ResFinder/all.faa CARD_filtered.faa ResFinderFG_unique.faa OS_unique.faa > ResFinder+CARD+FG+OS.faa
diamond makedb --in ResFinder+CARD+FG+OS.faa --db ResFinder+CARD+FG+OS --ignore-warnings
cp ResFinder+CARD+FG+OS.faa CLEVER.variants.faa
mkdir VARIANTS
mv *faa VARIANTS/
mv *dmnd VARIANTS/
mv *blastp VARIANTS/

## Match and cluster FAMILIES database (90% identity)
echo "Creating families database..."
#diamond makedb --in dbs/ResFinder/all.faa --db ResFinder --ignore-warnings
## If this is an update to the database, do this:
diamond makedb --in dbs/CLEVER/CLEVER.families.faa --db CLEVER --ignore-warnings
diamond blastp -q dbs/ResFinder/all.faa --db CLEVER -o ResFinder_vs_CLEVER.blastp --id 90 --query-cover 90 --subject-cover 90 --outfmt 6 --un ResFinder_unique.faa
cat dbs/CLEVER/CLEVER.families.faa ResFinder_unique.faa > ResFinder.faa
diamond makedb --in ResFinder.faa --db ResFinder --ignore-warnings
diamond blastp -q dbs/CARD/protein_fasta_protein_homolog_model.fasta --db ResFinder -o CARD_vs_ResFinder.blastp --id 90 --query-cover 90 --subject-cover 90 --outfmt 6 --un CARD_unique.faa
diamond makedb --in dbs/CARD/not_ARGs.fasta --db CARD_NOT_ARGs --ignore-warnings
diamond blastp -q CARD_unique.faa --db CARD_NOT_ARGs -o CARD_vs_CARD_NON_ARGs.blastp --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un CARD_filtered.faa
cd-hit -i CARD_filtered.faa -o CARD_clustered.faa -c 0.9 -n 5 -aS 0.9 -aL 0.9 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa > ResFinder+CARD.faa
diamond makedb --in ResFinder+CARD.faa --db ResFinder+CARD --ignore-warnings
diamond blastp -q dbs/ResFinderFG/ResFinder_FG_AA.faa --db ResFinder+CARD -o ResFinderFG_vs_ResFinder+CARD.blastp  --id 90 --query-cover 90 --subject-cover 90 --outfmt 6 --un ResFinderFG_unique.faa
cd-hit -i ResFinderFG_unique.faa -o ResFinderFG_clustered.faa -c 0.9 -n 5 -aS 0.9 -aL 0.9 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa ResFinderFG_clustered.faa > ResFinder+CARD+FG.faa
diamond makedb --in ResFinder+CARD+FG.faa --db ResFinder+CARD+FG --ignore-warnings
diamond blastp -q Other_sources.faa --db ResFinder+CARD+FG -o OS_vs_ResFinder+CARD+OS.blastp --id 90 --query-cover 90 --subject-cover 90 --outfmt 6 --un OS_unique.faa
cd-hit -i OS_unique.faa -o OS_clustered.faa -c 0.9 -n 5 -aS 0.9 -aL 0.9 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa ResFinderFG_clustered.faa OS_clustered.faa > ResFinder+CARD+FG+OS.faa
diamond makedb --in ResFinder+CARD+FG+OS.faa --db ResFinder+CARD+FG+OS --ignore-warnings
cp ResFinder+CARD+FG+OS.faa CLEVER.families.faa
mkdir FAMILIES
mv *faa FAMILIES/
mv *dmnd FAMILIES/
mv *blastp FAMILIES/
mv *clstr FAMILIES/

## Match and cluster LINEAGES database (70% identity)
echo "Creating lineages database..."
#diamond makedb --in dbs/ResFinder/all.faa --db ResFinder --ignore-warnings
## If this is an update to the database, do this:
diamond makedb --in dbs/CLEVER/CLEVER.lineages.faa --db CLEVER --ignore-warnings
diamond blastp -q dbs/ResFinder/all.faa --db CLEVER -o ResFinder_vs_CLEVER.blastp --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un ResFinder_unique.faa
cat dbs/CLEVER/CLEVER.lineages.faa ResFinder_unique.faa > ResFinder.faa
diamond makedb --in ResFinder.faa --db ResFinder --ignore-warnings
diamond blastp -q dbs/CARD/protein_fasta_protein_homolog_model.fasta --db ResFinder -o CARD_vs_ResFinder.blastp --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un CARD_unique.faa
diamond makedb --in dbs/CARD/not_ARGs.fasta --db CARD_NOT_ARGs --ignore-warnings
diamond blastp -q CARD_unique.faa --db CARD_NOT_ARGs -o CARD_vs_CARD_NON_ARGs.blastp --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un CARD_filtered.faa
cd-hit -i CARD_filtered.faa -o CARD_clustered.faa -c 0.7 -n 4 -aS 0.7 -aL 0.7 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa > ResFinder+CARD.faa
diamond makedb --in ResFinder+CARD.faa --db ResFinder+CARD --ignore-warnings
diamond blastp -q dbs/ResFinderFG/ResFinder_FG_AA.faa --db ResFinder+CARD -o ResFinderFG_vs_ResFinder+CARD.blastp  --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un ResFinderFG_unique.faa
cd-hit -i ResFinderFG_unique.faa -o ResFinderFG_clustered.faa -c 0.7 -n 4 -aS 0.7 -aL 0.7 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa ResFinderFG_clustered.faa > ResFinder+CARD+FG.faa
diamond makedb --in ResFinder+CARD+FG.faa --db ResFinder+CARD+FG --ignore-warnings
diamond blastp -q Other_sources.faa --db ResFinder+CARD+FG -o OS_vs_ResFinder+CARD+FG.blastp --id 70 --query-cover 70 --subject-cover 70 --outfmt 6 --un OS_unique.faa
cd-hit -i OS_unique.faa -o OS_clustered.faa -c 0.7 -n 4 -aS 0.7 -aL 0.7 -d 0
cat ResFinder_clustered.faa CARD_clustered.faa ResFinderFG_clustered.faa OS_clustered.faa > ResFinder+CARD+FG+OS.faa
diamond makedb --in ResFinder+CARD+FG+OS.faa --db ResFinder+CARD+FG+OS --ignore-warnings
cp ResFinder+CARD+FG+OS.faa CLEVER.lineages.faa
mkdir LINEAGES
mv *faa LINEAGES/
mv *dmnd LINEAGES/
mv *blastp LINEAGES/
mv *clstr LINEAGES/


## Collect databases
echo "Collecting databases..."
cp VARIANTS/CLEVER.variants.faa .
cp FAMILIES/CLEVER.families.faa .
cp LINEAGES/CLEVER.lineages.faa .

## Check which genes exist on plasmids
echo "Building plasmid database..."
makeblastdb -dbtype nucl -out PLSDB -in dbs/PLSDB/plsdb.fna
echo "Identifying mobile ARG variants..."
tblastn -db PLSDB -query CLEVER.variants.faa -out CLEVER.variants.vs.PLSDB.tblastn -outfmt 6 -num_threads 64 -evalue 1e-30
echo "Identifying mobile ARG families..."
tblastn -db PLSDB -query CLEVER.families.faa -out CLEVER.families.vs.PLSDB.tblastn -outfmt 6 -num_threads 64 -evalue 1e-30
echo "Identifying mobile ARG lineages..."
tblastn -db PLSDB -query CLEVER.lineages.faa -out CLEVER.lineages.vs.PLSDB.tblastn -outfmt 6 -num_threads 64 -evalue 1e-30

echo "Filtering results..."
perl build-scripts/filter_blast.pl CLEVER.lineages.vs.PLSDB.tblastn 70
perl build-scripts/filter_blast.pl CLEVER.families.vs.PLSDB.tblastn 90
perl build-scripts/filter_blast.pl CLEVER.variants.vs.PLSDB.tblastn 98

## Make gene annotations
echo "Annotating ARGs..."
diamond makedb --in dbs/CLEVER/CLEVER.families.faa --db CLEVER.families --ignore-warnings
diamond blastp -q CLEVER.variants.faa --db CLEVER.families -o CLEVER_vs_families.blastp --id 90 --query-cover 90 --subject-cover 90 --outfmt 6
perl build-scrpts/annotate_args.pl CLEVER_vs_families.blastp blacklist.txt $VERSION

perl build-scrpts/rename_fasta.pl CLEVER.variants.faa $CLEVER_ANNOTATION CLEVER.variants.vs.PLSDB.tblastn.hits.txt $VERSION > CLEVER.variants.final.faa
perl build-scrpts/rename_fasta.pl CLEVER.families.faa $CLEVER_ANNOTATION CLEVER.families.vs.PLSDB.tblastn.hits.txt $VERSION > CLEVER.families.final.faa
perl build-scrpts/rename_fasta.pl CLEVER.lineages.faa $CLEVER_ANNOTATION CLEVER.lineages.vs.PLSDB.tblastn.hits.txt $VERSION > CLEVER.lineages.final.faa

## Finalize CLEVER build
echo "Finalizing CLEVER build..."
cat CLEVER.variants.final.faa | grep ">" | sed "s/^>//" | sed "s/ /\t/g" > CLEVER.variants.tsv
cat CLEVER.families.final.faa | grep ">" | sed "s/^>//" | sed "s/ /\t/g" > CLEVER.families.tsv
cat CLEVER.lineages.final.faa | grep ">" | sed "s/^>//" | sed "s/ /\t/g" > CLEVER.lineages.tsv

mv CLEVER.variants.final.faa CLEVER-build/CLEVER.variants.faa
mv CLEVER.families.final.faa CLEVER-build/CLEVER.families.faa
mv CLEVER.lineages.final.faa CLEVER-build/CLEVER.lineages.faa
mv CLEVER.*.tsv CLEVER-build/
mv CLEVER_vs_families.blastp.mapping.txt CLEVER-build/CLEVER.variant-family-mapping.txt

tar -czvf CLEVER.tgz CLEVER-build/

echo "Done!"
echo "Your new build of CLEVER is located in CLEVER-build/"
