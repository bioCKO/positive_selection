#XP_015777750.1
#PREDICTED: GDP-mannose 4,6 dehydratase-like isoform X2
#associated with GO term:
	go='GO:0007219' #notch signaling pathway


N sepcies = 5
alignment length = 286 amino acids

From Uniprot:
	Function:
	Catalyzes the conversion of GDP-D-mannose to GDP-4-dehydro-6-deoxy-D-mannose.

	The connection between GDP-mannose 4,6, dehydratase and notch seems weak.
	Uniprot association is 'Inferred from Sequence or Structural Similarity'
	One article that seems to address it is:
	Smith et et. 2002: Conditional control of selectin ligand expression and global fucosylation events in mice with a targeted mutation at the FX locus
	This article was cited in:
	Ishikawa et al. 2010: Two Pathways for Importing GDP-fucose into the Endoplasmic Reticulum Lumen Function Redundantly in the O-Fucosylation of Notch in Drosophila*



From Hemmond et al. 2014:
"In our dataset, notch1 and notch2 transmembrane protein genes 
and a regulatory gene, GDP-mannose 4,6, dehydratase, were consistently 
up-regulated in tips, while one Notch regulatory gene, E3 
ubiquitin-protein ligase MIB2, was up-regulated in bases."







Naive Empirical Bayes (NEB) analysis (please use the BEB results.)
#!note did not use the BEB results here because nearly all sites were listed, but all with posterior probabilities ~0.5
Positive sites for foreground lineages Prob(w>1):

     1 K 0.506
     3 I 0.945
    25 D 0.571
    52 L 0.978*  V>L
    57 I 0.988*  L>I
    68 H 0.918
    69 R 0.610
    70 I 0.509
    85 R 0.882
    92 N 0.896
   120 K 0.771
   167 A 0.997**  S>A
   189 N 0.676
   193 Y 0.996**  F>Y
   207 T 0.971*   V>T
   226 S 0.943
   234 N 0.991**  T>N
   236 N 0.770
   240 A 0.834
   266 N 0.992**  G>N
   268 V 0.994**  K>V
   270 E 0.859
   285 V 0.534



