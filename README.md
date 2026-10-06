PHYLOGENOMIC ANALYSIS OF DIPTERA RELATIONSHIPS: RESOLVING CONFLICT IN RAPID RADIATION EVENTS

Project Overview
================

This repository contains the complete bioinformatics workflow for phylogenomic reconstruction 
of relationships among major Diptera clades (Drosophilidae, Calyptratae, and other families). 
The analysis uses 178 single-copy orthologous genes across 137 species to investigate 
evolutionary relationships and resolve long-standing phylogenetic conflict caused by rapid 
radiation events and incomplete lineage sorting (ILS).

This work directly builds upon foundational studies (Wiegmann et al., 2011) while addressing 
limitations in prior phylogenetic inference through comprehensive genomic sampling and 
advanced analytical approaches.

Full thesis available: https://hdl.handle.net/20.500.14239/31786


Key Objectives
==============

1. Reassess phylogenetic relationships within Diptera using phylogenomic approaches
2. Evaluate effects of rapid radiation events and incomplete lineage sorting on inference
3. Compare concatenation vs coalescent-based phylogenetic methods
4. Assess robustness of inferred relationships through multiple analytical strategies
5. Optimize computational workflows for large-scale phylogenomic analyses on HPC


Dataset & Methods
=================

Species & Genes:
- 137 species across Diptera
- 178 single-copy orthologous genes (BUSCO protein set)
- Multiple sequence alignments per gene

Analytical Approaches:
- Concatenation-based inference (IQ-TREE, RAxML-NG)
- Coalescent-based inference (ASTRAL)
- Amino-acid, nucleotide, and codon-based datasets
- Bootstrap support evaluation (1000 replicates)
- Model selection testing


COMPLETE WORKFLOW & COMMANDS
=============================

STEP 1: BUSCO GENE EXTRACTION & PREPARATION
--------------------------------------------

Identify and extract conserved BUSCO genes from genomic assemblies

**Command line:**
```bash
busco -i [species.genome.fasta] -l diptera_odb12 \
-o [species_buscos] -m genome -c 6 --metaeuk
```

Output: Single copy BUSCO genes for each species   
Result: 178 single-copy BUSCO genes present in all the 137 species


STEP 2: ALIGNMENT AND TRIMMING
--------------------------------------------

Align each of the 178 orthologous BUSCO gene sets  
Perform multiple sequence alignment for each gene separately to remove poorly aligned regions  

Tools used:   
MAFFT v7.526  
trimAl v1.5.rev0  
   
Output: 178 aligned gene sequences (one file per gene)


STEP 3A: CONCATENATION-BASED PHYLOGENETIC INFERENCE (IQ-TREE)
--------------------------------------------------------------

Concatenate all 178 genes, gene-by-gene for each species to reconstruct a single supermatrix  
Reduced datasets were also created with less genes, less species, and a subset of "good genes"

For all datasets: 

**Run IQ-TREE with  model selection:**
```bash
iqtree -s [alignment.fasta] -m [substitution_model] -B 1000 -T AUTO -pre [output_ prefix_name]
```

Parameters:
-s: specifies input sequence alignment file in FASTA format   
-m MFP: Model finder  (tests multiple models)    
-B 1000: 1000 ultrafast bootstrap replicates   
-T AUTO: Automatically use available CPU threads   

Output: 
- .treefile: Best phylogenetic tree
- .log: Analysis log with model selection results
- .contree: Consensus tree with bootstrap support


STEP 3B: CONCATENATION-BASED INFERENCE (RAxML-NG)
---------------------------------------------------

Alternative concatenation analysis for robustness testing

**Command line:**
```bash
raxml-ng --all --msa [alignment.fasta] --model [LG+F+I+G] --prefix [Diptera_phylogeny] --seed 693044 --bs-metric tbe --tree rand{1} --bs-trees 100
```

Substitution model for amino acid dataset: LG+F+I+G   
Substitution model for nucleotide dataset: GTR+I+G  

Parameters:  
--msa: Species multiple sequence alignment FASTA file   
--model: Substiution model  
--seed 693044: Random seed selected to ensure reproducibility 
--bs-trees 100: 100 bootstrap replicates  

Output: 
- .bestTree: Best-scoring maximum likelihood tree
- .support: Bootstrap support values


STEP 4: GENE TREE INFERENCE (INDIVIDUAL GENES)
-----------------------------------------------

Build phylogenetic tree for each individual gene (required for ASTRAL)  

For each of 178 genes:  
iqtree command line was ran  

**Collect all gene trees:**
```bash
cat gene_*.treefile > all_gene_trees.txt
```

Then analysis was conducted using weighted ASTRAL (wASTRAL), implemented in the package ASTER  

Output: 178 individual gene tree files (.treefile format)


STEP 5: COALESCENT-BASED SPECIES TREE INFERENCE (ASTRAL)
---------------------------------------------------------

Infer species tree from gene trees using ASTRAL algorithm

Output: 
- astral_species_tree.tre: Species tree with local posterior probabilities
- ASTRAL handles ILS and gene tree discordance


STEP 6: TREE COMPARISON & TOPOLOGY ANALYSIS
---------------------------------------------

Compare topologies from different inference methods

**Calculate Robinson-Foulds distances between trees:**
```bash
Rscript compare_trees.R iqtree_protein.treefile raxml_protein.bestTree astral_species_tree.tre
```

Compare support values across methods:  
- Bootstrap support (IQ-TREE, RAxML-NG): ≥95% = strong support  
- Local posterior probabilities (ASTRAL): ≥0.95 = strong support  
- Identify conflicting nodes between concatenation and coalescent approaches 


STEP 7: HPC OPTIMIZATION & RUNTIME REDUCTION
----------------------------------------------

Initial runs (local machine / small cluster):
- 178 gene alignments: ~2-3 weeks
- IQ-TREE concatenated analysis: ~5-7 days
- Full RAxML-NG analysis: ~4-6 days

Optimizations for HPC:
Parallelize gene tree inference across nodes:
   - Submit 178 gene tree jobs as array jobs
   - Each gene processes independently
   - Reduces 2-3 weeks to 1-2 days


Request appropriate resources:
   - 16+ CPU cores
   - 64GB+ RAM
   - High-speed storage for temporary files


KEY RESULTS
===========

Phylogenetic Topology:
- Resolved relationships among major Diptera clades with a ((Drosophilidae, Calyptratae), Tephritidae)) configuration
- Evaluated support for competing hypotheses
- Identified nodes with strong vs weak support across methods

Incomplete Lineage Sorting:
- Quantified discordance between gene trees
- ASTRAL posterior probabilities indicate degree of ILS
- Short internal branches support hypothesis of rapid radiation

Dataset Comparison:
- Amino-acid vs nucleotide vs codon alignments show similar topologies
- Concatenation vs coalescent approaches generally congruent
- Identified regions of conflict related to ILS

Biological Interpretation:
- Results inform understanding of Diptera evolutionary history
- Support for rapid diversification events in early Diptera evolution
- Future work: increased taxon sampling, modeling gene flow


REPRODUCIBILITY & DEPENDENCIES
===============================

Software Requirements:
- MAFFT v7.526 (sequence alignment)
- trimAl v1.5.rev0 (alignment trimming)
- IQ-TREE v3.0.1 (phylogenetic inference)
- RAxML-NG v1.2.2 (alternative phylogenetic inference)
- ASTRAL (in the package ASTER) (coalescent-based species tree)
- figtree v1.4.4 (tree visualization)
- iTOL v7 (tree visualisation)
- R (tree comparision)

Computing Environment:
- HPC cluster with SLURM job scheduler
- Minimum 16 CPU cores, 64GB RAM


Version Control:
All scripts and workflows version-controlled for reproducibility


HOW TO REPRODUCE THIS ANALYSIS
===============================

1. Obtain genomic data for  species
2. Run BUSCO to extract single-copy orthologous genes
3. Execute MAFFT alignment for each gene
4. Trim alignments with trimAl
5. Run both concatenation methods (IQ-TREE & RAxML-NG)
6. Infer individual gene trees (IQ-TREE & RaxML-NG)
7. Run ASTRAL analysis 
8. Compare results
9. Optimize on HPC for faster runtime

Expected output: Robust phylogenetic hypothesis for Diptera relationships


BIOLOGICAL SIGNIFICANCE
=======================

This study addresses:
- Decades-long conflict in Diptera phylogeny (since Wiegmann et al., 2011)
- Challenges of phylogenetic inference under rapid radiation
- Role of incomplete lineage sorting in deep divergences
- Advantages of phylogenomic approaches for resolving ancient radiations

Findings contribute to:
- Better understanding of Diptera evolution
- Methodological insights for phylogenomics in rapid radiation contexts
- Framework for future studies on Diptera and other rapidly diversifying groups


REFERENCES & CITATIONS
======================

Wiegmann BM, et al. (2011). Episodic radiations in the fly tree of life. 

Full thesis with detailed literature review:
https://hdl.handle.net/20.500.14239/31786

For detailed methods: See thesis chapters 3.2-3.8


AUTHOR & ACKNOWLEDGMENTS
========================

Analysis conducted by: Vani Jain
Institution: Università di Pavia
Date: February - November 2025
Thesis Grade: 30/30

