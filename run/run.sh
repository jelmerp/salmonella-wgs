#!/bin/bash

# ==============================================================================
# ASSEMBLY & ANNOTATION
# ==============================================================================
## Run Trimmomatic
for R1 in data/fastq/*_R1_*fastq.gz; do
    sbatch mcic-scripts/trim/trimmomatic.sh \
        -i "$R1" -o results/trim \
        -a mcic-scripts/trim/adapters.fa \
        -p "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" \
        -A "2:30:10"
done

## Run SPAdes - #? kmer-sizes from Carroll et al 2021
for R1 in results/trim/*_R1_*fastq.gz; do
    sampleID=$(basename "$R1" | sed 's/_.*//')
    sbatch mcic-scripts/assembly/spades.sh \
        -i "$R1" -o results/spades/"$sampleID" \
        -m isolate -k "21,33,55,77,99,127"
done

## Quickly check some assembly stats
conda activate /users/PAS0471/jelmer/miniconda3/envs/bbmap-env
shopt -s globstar
for assembly in results/spades/*/contigs.fasta; do
    stats.sh "$assembly"
done

## Run BBmap to estimate the mean coverage for each sample
for R1 in results/trim/*_R1_*.fastq.gz; do
    R2=${R1/_R1_/_R2_}
    sampleID=$(basename "$R1" | sed 's/_S.*//')
    assembly=results/spades/"$sampleID"/contigs.fasta
    outfile=results/bbmap/"$sampleID".bam
    sbatch scripts/bbmap.sh "$R1" "$R2" "$assembly" "$outfile"
done

## Copy & rename the assembly FASTAs for easier access
mkdir -p results/spades/all_assemblies
for dir in results/spades/AG*; do
    cp -v "$dir"/contigs.fasta results/spades/all_assemblies/"$(basename "$dir")".fasta
done

## Run QUAST to check genome quality
for assembly in results/spades/all_assemblies/*fasta; do
    sbatch mcic-scripts/assembly/quast.sh -i "$assembly" -o results/quast
done

## Run checkM to check genome quality -- for all assemblies at once
sbatch mcic-scripts/assembly/checkm.sh -i results/spades/all_assemblies -o results/checkm

## Run Kraken to check for contamination
kraken_db=/fs/project/PAS0471/jelmer/refdata/kraken/std-plus-fungi
for asm in results/spades/all_assemblies/*fasta; do
    sbatch mcic-scripts/metagenomics/kraken-run.sh -i "$asm" -o results/kraken -d "$kraken_db" -n
done

## Annotate the assemblies with Prokka
for assembly in results/spades/all_assemblies/*fasta; do
    sbatch mcic-scripts/annot/prokka.sh -i "$assembly" -o results/prokka -g Salmonella -s enterica
done

# ==============================================================================
# PANGENOMICS
# ==============================================================================
## Roary
sbatch mcic-scripts/bact/roary.sh -i tmp_gff -o results/roary


# ==============================================================================
# PLASMIDS & PHAGES
# ==============================================================================
## Look for plasmids with PlasmidFinder
for asm in results/spades/all_assemblies/*fasta; do
    outdir=results/plasmidfinder/$(basename "$asm" .fasta)
    sbatch mcic-scripts/misc/plasmidfinder.sh -i "$asm" -o "$outdir"
done

## Phaster #! Doesn't work - jobs in queue at server for days
software/phaster_scripts/phaster.py --fasta results/spades/all_assemblies/AG21-0050.fasta
software/phaster_scripts/phaster.py --get-status


# ==============================================================================
# ANTIMICROBIAL RESISTANCE & VIRULENCE
# ==============================================================================
## Check for antimicrobial resistance with Amrfinderplus
for assembly in results/spades/all_assemblies/*fasta; do
    aa_fasta=results/prokka/$(basename "$assembly" .fasta).faa
    gff=results/prokka/$(basename "$assembly" .fasta).gff
    sbatch mcic-scripts/misc/amrfinderplus.sh -n "$assembly" -p "$aa_fasta" -g "$gff" -o results/amrfinderplus -s Salmonella
done

## Check for virulence with Abricate
# source activate /fs/project/PAS0471/jelmer/conda/abricate-1.0.1


# ==============================================================================
# PHYLOGENETICS
# ==============================================================================
## Get ref assembly FASTA and GFF from accession numbers
sbatch scripts/download-refs.sh results/patric_refs/input/ref_IDs.txt results/patric_refs

## Run kSNP3 to get an alignment of core SNPs -- own genomes only
mkdir -p results/ksnp3/own_asm
paste <(ls "$PWD"/results/spades/all_assemblies/*) <(ls results/spades/all_assemblies/ | sed 's/.fasta//') > results/ksnp3/own_asm/assembly_list.txt
sbatch mcic-scripts/trees/ksnp3.sh -i results/ksnp3/own_asm/assembly_list.txt -o results/ksnp3/own_asm

## Run kSNP3 to get an alignment of core SNPs -- own genomes and selected refs
mkdir -p results/ksnp3/with_refs/refs
for fna in results/patric_refs/*fna; do
    ## Rename FASTA files for knsp3 ... no extra periods, etc
    newname=results/ksnp3/with_refs/refs/$(basename "$fna" | sed 's/\..*//').fasta
    cp -fv "$fna" "$newname"
done
cp results/ksnp3/own_asm/assembly_list.txt results/ksnp3/with_refs/assembly_list.txt
paste <(ls "$PWD"/results/ksnp3/with_refs/refs/*fasta) <(ls results/ksnp3/with_refs/refs | sed 's/.fasta//') >> results/ksnp3/with_refs/assembly_list.txt
sbatch mcic-scripts/trees/ksnp3.sh -i results/ksnp3/with_refs/assembly_list.txt -o results/ksnp3/with_refs

## Plot the resulting tree
Rscript mcic-scripts/trees/ggtree.R -i results/ksnp3/with_refs/tree.core.tre -o results/ksnp3/with_refs/tree.core.png

## Run IQ-tree to create a tree based on the kSNP3 SNPs #! fails
sbatch --mem=4G mcic-scripts/trees/iqtree.sh -i results/ksnp3/core_SNPs_matrix.fasta -p results/iqtree/knsp3/noboot/all_core
sbatch -t 24:00:00 mcic-scripts/trees/iqtree.sh -i results/ksnp3/core_SNPs_matrix.fasta -p results/iqtree/knsp3/boot1k/all_core -b 1000 
#! ERROR: phylotree.cpp:2730: virtual double PhyloTree::optimizeAllBranches(int, double, int): Assertion `fabs(new_tree_lh-tree_lh) < max_delta_lh' failed.


# ==============================================================================
#                          TIMETREES - CERRO ONLY
# ==============================================================================
clock_rate=2.4e-7  #2.4E-7/site/year (from Rodriguez-Rivera et al., 2014, Carole et al)
sero=cerro
ref_id="AG21-0559"
dates_master=data/meta/dates.tsv
dates=results/dates/"$sero"_dates.tsv && mkdir -p results/dates

## Prepare dates file
awk -v sero=$sero 'NR==1 || tolower($5) == sero' "$dates_master" | cut -f 1,9 |
    sed "s/$ref_id/Reference/" > "$dates" # Snippy will rename ref. seq. as "Reference"

## Core SNPs with Snippy => Align FASTQs to ref
genome_list=results/serotype_lists/"$sero".txt
ref=results/spades/all_assemblies/"$ref_id".fasta
outdir=results/snippy/"$sero" && mkdir -p "$outdir"
tab="$outdir"/"$sero"_input.tsv
find data/fastq | sort | grep -f <(sort "$genome_list") |
    paste <(sort "$genome_list") - - | grep -v "$ref_id" > "$tab"
sbatch mcic-scripts/bact/snippy-multi.sh -i "$tab" -r "$ref" -o "$outdir"

## Remove recombinant regions with Gubbins - note that it requires whole-genome alignments/invariant sites
aln=results/snippy/"$sero"/core.full.aln
sbatch mcic-scripts/bact/gubbins.sh -i "$aln" -o results/gubbins/"$sero"/"$sero"

## Try Gubbins with dates
sbatch mcic-scripts/bact/gubbins.sh -i "$aln" -o results/gubbins/"$sero"/"$sero" -d "$dates"

## Treetime
fa=results/gubbins/"$sero"/"$sero".filtered_polymorphic_sites.fasta
tree=results/gubbins/"$sero"/"$sero".final_tree.tre
### A. Without fixed clock rate
sbatch mcic-scripts/trees/treetime.sh -i "$fa" -d "$dates" -t "$tree" -o results/treetime/"$sero"
cat results/treetime/"$sero"/dates.tsv
### B. With fixed clock rate
sbatch mcic-scripts/trees/treetime.sh -i "$fa" -d "$dates" -t "$tree" -c "$clock_rate" -o results/treetime/"$sero"_fixedrate
cat results/treetime/"$sero"_fixedrate/dates.tsv

## Detect Temporal signal with TempEst - followed https://beast.community/tempest_tutorial
#> Date range	             2.2809
#> Slope (rate)	             -105.4153
#> X-Intercept (TMRCA)	     2021.6657
#> Correlation Coefficient	 -0.9887
#> R squared	             0.9774
#> Residual Mean Squared	 109.2983
module load R/4.2.1-gnu11.2
R --quiet -e "1 - pt(0.9887 * sqrt( (14 / (1 - 0.9887^2) ) ), df = 14)" #> 3.061995e-13
#! Note that I got a different result using the BactDating R package, see scripts/bactdating.R
