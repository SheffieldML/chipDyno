CHIPDYNO software
Version 0.1.0		Friday 14 Mar 2008 at 17:11

This code contains the scripts for recreating the results in G.Sanguinetti, M.
Rattray and Neil D. Lawrence, "A probabilistic dynamical model for quantitative inference of the regulatory mechanism of transcription".

Version 0.121
-------------
Add some files that were missing from the previous release.

Version 0.12
------------
Version 0.11
------------
Release of code for final version of Bioinformatics paper.
Version 0.1
-----------
Initial release.

R Files
------------
R scripts associated with the toolbox are: (Sorted)
chipChipTextRead.R: reads TXT file for the Lee ChIP data files.
chipDynoActTransFact.R: identifies significantly varying TFs.
chipDynoActTransFactNoise.R: identifies significantly varying TFs.
chipDynoCompareTrans.R: compares list of active tfs at different p
chipDynoCompareX.R: compares connectivity for different pvals
chipDynoExpectationsFast.R: computes posterior expectations of TFA.
chipDynoExpectationsFastNoise.R: computes posterior expectations TFA.
chipDynoGeneAct.R: given a gene, lists activators in decreasing order
chipDynoGeneActNoise.R: given a gene, lists activators in decreasing order
chipDynoLikeStat.R: marginal likelihood for chipChip dynamical model
chipDynoLikeStatGrad.R: gradient of CHIPDYNOLIKESTAT.
chipDynoLikeStatNoise.R: marginal likelihood for chipChip dynamical model
chipDynoLikeStatNoiseGrad.R: gradient of CHIPDYNOLIKESTATNOISE.
chipDynoLoadData.R: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap1.R: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap1p2.R: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap2.R: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap2p5.R: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap5.R: loads Spellman Data with Lee et al ChIP data.
chipDynoMaxDiff.R: computes most significant changes in TFAs
chipDynoMultRegs.R: finds genes that have multiple regulators
chipDynoNewTransFactNoise.R: tfs active for us and not for Tu et al.
chipDynoPerGenes.R: genes periodic according to Tu et al
chipDynoSignRels.R: finds the regulatory strength conserved rels
chipDynoStatPostEst.R: computes posterior expectations 
chipDynoStatPostEstNoise.R: computes posterior expectations.
chipdynoToolboxes.R: The toolboxes needed for CHIPDYNO software.
chipDynoTransFact.R: provides gene-specific TFAs with errorbars.
chipDynoTransFactCorrel.R: correlations between transcription factors
chipDynoTransFactNoise.R: given a transcription factor, provides TFAs.
chipDynoTuLoadData.R: loads Tu Data with Lee et al ChIP data.
chipReduceVariables.R: reduce  number of variables in chipDyno model
chipStatMatrixInverter.R: inverts block tridiagonal matrices for chipChip
chipStatMatrixInverterNoise.R: inverts block tridiagonal matrices for chipChip
chipTextRead.R: reads TXT file for the Spellman data files.
chipTuTextRead.R: assigns common names to probe IDs
demSpellmanDynoStat.R: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanp1.R: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanp1p2.R: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanp2.R: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanp5.R: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanp5p2.R: demonstrates dynamical chipCHIP on Spellman data.
demTu.R: demonstrates dynamical chipCHIP on Tu data.
demTuNoNoise.R: demonstrates dynamical chipCHIP on Tu data.
