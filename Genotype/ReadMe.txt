This directory contains the genotype data for the mice involved in a Jackson Laboratory Aging 
Center (JAC) study that ran from 2011 to 2016. The study consisted of 600 (300 female and 300 
male) Diversity Outbred (DO) mice in a longitudinal cohort and another 543 (271 female and 272 
male) DO mice split over three cross-sectional cohorts that were sacrificed for tissue collection 
at approximately 6, 12, or 18 months. They were genotyped using a MegaMUGA array.


The folder contains 5 files:
	-	QTL2_cross2_LongAndCS.rds: The Q/qtl2::cross2 data object after the cleaning/processing described 
		above. It contains 1097 samples/mice (from both the longitudinal and cross-sectional study arms)
		and 70013 markers across the 19 autosomes, the sex chromosomes (2143 markers on the X and 35 on
		the Y), the mitochondrial DNA (32 markers), and the pseudo-autosomal region (2 markers)
	-	Long_GenoProbs.rds: The 36-state genoprobs list (created using R/qtl2::calc_genoprobs) for mice
		in the longitudinal arm
	-	Long_AlleleProbs.rds: The 8-state alleleprobs list (created using R/qtl2::genoprob_to_alleleprob)
		for mice in the longitudinal arm
	-	CS_GenoProbs.rds: The 36-state genoprobs list (created using R/qtl2::calc_genoprobs) for mice
		in the cross-sectional arm
	-	CS_AlleleProbs.rds: The 8-state alleleprobs list (created using R/qtl2::genoprob_to_alleleprob)
		for mice in the cross-sectional arm


These genotypes have been processed, during which the following changes were made:
	Regarding samples (mice), 4 sample mix-ups (2 from the longitudinal arm and 2 from the 
	cross-sectional arm) that were identified by coat color mismatches were fixed and 16 samples 
	were dropped:
		-	9 mice (3 longitudinal, 6 cross-sectional) that were missing 20% of their genotype data 
			OR had a genotyping error rate of 0.3% OR were missing 10% and had an error rate of 0.25%. 
		- 	6 mice (4 longitudinal, 2 cross-sectional) that were involved in duplications (identical 
			genotype data to another mouse)
		- 	4 mice (2 longitudinal, 2 cross-sectional) with irreconcilable coat color mismatches (recorded 
			coat color does not match coat color genotype)

	Regarding markers, a total of 7712 (~10%) out of 77725 markers were dropped:	
		- 	737 markers that were missing genotype data for every sample
		- 	5336 markers for which the founders have identical genotypes and are thus non-informative
		- 	1639 markers that were missing for 20% of the samples and/or they had a genotyping error 
			rate of 5%


The raw data and data processing files are available on tier2 data storage in:
	churchilldev.jax.org::/tier2/churchill-lab/JAC/data/genotypes_agingDO_long_and_cs/



