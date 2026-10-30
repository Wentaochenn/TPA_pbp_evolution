# Figure legends

## Figure 1. Global distribution and PBP variation in the Treponema pallidum collection

A, Geographical distribution of 1,861 genomes from 31 countries; symbol area represents genome count. B–D, Country-specific genome counts, collection-year distributions and Nichols–SS14 proportions, respectively. Countries follow a common order across B–D; symbol area in C represents genome count. E, Carrier counts for variants in *mrcA*, *pbp1*, *pbp2* and *tp47*, with *lytB* included as a pathway-related comparator. Labels give carriers divided by genomes with a callable genotype for that variant. F, Collection-year distributions of variant carriers; symbol area represents carrier count. Crimson symbols identify the six focal PBP variants and grey symbols the remaining displayed variants. Only genomes with callable genotypes and available collection years contribute to F. Variant-specific denominators exclude missing calls.

## Figure 2. Lineage distribution and recurrent emergence of PBP variants

A, Six focal substitutions mapped onto the global phylogeny of 1,861 genomes. Blue and crimson sectors indicate Nichols and SS14, respectively. Coloured tips identify carriers; coloured internal branches connect descendant clades composed entirely of observed carriers. B, Country-level variant carriage versus the proportion of SS14 among the same variant-specific callable genomes. Countries with at least three callable genomes are labelled; symbol area represents callable genome count. Solid lines show unweighted linear regressions, r denotes Pearson correlation, and dashed diagonals indicate equality of the percentages. The upper axis gives the complementary Nichols proportion. C, Lineage-specific carriage percentages with carrier/callable counts. D, Minimum and maximum numbers of gains across equally parsimonious, equal-cost histories with an unrestricted root state. E, Gain and loss counts from 100 stochastic maps under the equal-rates model, using known-state tips for each variant. Orange circles and green triangles indicate gains and losses; points and intervals represent means and the 2.5th–97.5th percentiles. Parsimony bounds in D describe equally optimal histories, whereas intervals in E describe the mapping distributions.

## Figure 3. Codon-level selection signals in PBP-associated genes

A, FUBAR posterior probabilities that the nonsynonymous substitution rate exceeds the synonymous rate at all 2,594 codons in *mrcA*, *pbp1*, *pbp2* and *tp47*. The dashed line marks a posterior probability of 0.90. The five focal codons are labelled and highlighted. The four gene alignments contain 1,861, 1,861, 1,860 and 1,859 genomes, respectively.

## Figure 4. MrcA protein types across the Nichols and SS14 lineages

A, Amino-acid states at variable positions in 16 MrcA protein types. WT denotes the wild-type sequence, and Mutation #1–#15 distinguish the remaining protein types. Grey cells retain the wild-type residue; purple cells carry an alternative residue. Protein types are resolved in 1,328 genomes, comprising 387 Nichols and 941 SS14 genomes. B, Protein-type counts by lineage. Symbol size represents genome count, with exact counts labelled; singleton types retain a minimum visible symbol size. C, The global phylogeny aligned with lineage and protein-type assignments. Blue-grey and purple marks indicate WT and variant-type membership, respectively. Each typed genome belongs to one protein type; the 533 untyped genomes remain on the tree and are identified separately. Numbers above protein-type columns indicate genome counts. The scale bar gives substitutions per site.

## Figure 5. Time-calibrated SS14 phylogenies and PBP acquisition histories

A,B, Maximum-clade-credibility phylogenies with common-ancestor node heights for the sampled 533-genome SS14 dataset under strict and uncorrelated lognormal relaxed clocks, respectively. Horizontal branch lengths represent calendar time. Terminal tracks show the six focal substitutions: coloured cells indicate carriers, light cells noncarriers and grey cells uncallable genotypes. Dashed vertical lines mark 1943. C, Dates of nine reciprocally matched gain events with matching support of at least 0.50 under both clocks. Purple circles and blue-grey triangles represent strict- and relaxed-clock medians; horizontal intervals show 95% highest posterior density intervals. Grey bubbles show carrier collection years for the corresponding strict-clock event, with area proportional to carrier count. Event numbers identify matched reconstructions within each substitution. D, Posterior probabilities that the matched acquisition dates occurred after 1943, shown separately for each clock; labels are rounded to three decimals.

# Supplementary figure legends

## Supplementary Figure 1. Transition-model sensitivity of reconstructed gain counts

Mean gain counts and 2.5th–97.5th percentiles from 100 stochastic maps under the equal-rates model (ER, blue) and all-rates-different model (ARD, red) for 23 PBP variants. ER constrains gain and loss rates to equality, whereas ARD estimates them separately. Points and horizontal segments show the means and percentile intervals, with numerical summaries labelled.

## Supplementary Figure 2. Structural environment of focal MrcA residues

A, AlphaFold reference structure of MrcA/TP0705 corresponding to UniProt O83703, rendered from the archived PDB coordinates using fixed PyMOL camera parameters. A506 is highlighted in red and M625 and G708 in orange. B, Enlarged view of their local structural environment. C, Per-residue predicted local distance difference test scores (pLDDT), extracted from the AlphaFold coordinate file; A506, M625 and G708 are labelled. D, Solvent-accessible surface areas (SASA) of MrcA residues A506, M625 and G708, calculated using 2,880 sphere points per atom and a 1.4 Å probe. Labels give residue-level areas in Å².

## Supplementary Figure 3. ESM-2 scores for PBP amino-acid variation

A, Position-level ESM-2 profiles for MrcA, Pbp1, Pbp2 and Tp47. Each position represents the arithmetic mean of log-likelihood ratios across its 19 non-reference amino acids. Open circles locate positions with observed substitutions, which are labelled; labels sharing a codon are combined. B, Substitution-specific scores for MrcA A506T, A506V, A506M, M625V and G708S, Pbp1 P564L and Tp47 S394R. The log-likelihood ratio is log P(alternative residue | masked sequence context) − log P(reference residue | masked sequence context). Negative and positive values indicate lower and higher conditional sequence probabilities than the reference residue, respectively. Scores are calculated in the Nichols reference-sequence context. Scores for all 40 observed PBP substitutions are reported in Supplementary Table 2.

## Supplementary Figure 4. Temporal signal across the genomic sampling schemes

Root-to-tip genetic distance versus collection year for the complete dated collection, the Nichols and SS14 lineages, and the candidate SS14 sampling subsets. Points represent genomes; lines show linear regressions and shaded bands their 95% confidence intervals. Blue indicates a positive slope and observed R² above the 95th percentile from date randomization; red indicates that at least one of these criteria is not met. Facet labels give the number of genomes, observed R², regression slope and randomized 95th-percentile R². These sample sizes identify the diagnostic datasets; the final BEAST2 dataset contains 533 SS14 genomes.

# Supplementary tables

## Supplementary Table 1. Genome metadata, amino-acid states and *mrcA* protein types

Sequence accessions, geographical and temporal annotations, lineage assignments, BioProject and BioSample identifiers, source studies, amino-acid states at 42 surveyed positions in *mrcA*, *pbp1*, *tp47*, *lytB* and *pbp2*, and *mrcA* protein types. Each position reports the observed residue or uncallable; the workbook footnotes specify the callability criteria. Protein types that could not be assigned are labelled untyped. Carrier and callable-genome totals are provided for each position.

## Supplementary Table 2. PBP amino-acid variants and substitution-specific ESM-2 scores

The 40 observed substitutions in *mrcA*, *pbp1*, *pbp2* and *tp47*, with locus identifiers, callable denominators, carriage frequencies in the full collection and each lineage, and substitution-specific ESM-2 scores.

## Supplementary Table 3. Sensitivity of focal FUBAR results to genome-wide phylogenetic topology

Primary-tree FUBAR posterior probabilities, their minimum, median and maximum across 20 sampled genome-wide bootstrap topologies, and the number of topologies with posterior probabilities exceeding 0.90.
