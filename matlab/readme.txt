CHIPDYNO software
Version 0.121		Friday 14 Mar 2008 at 17:11

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


MATLAB Files
------------

Matlab files associated with the toolbox are:

demSpellmanp1.m: demonstrates dynamical chipCHIP on Spellman data.
chipdynoToolboxes.m: The toolboxes needed for CHIPDYNO software.
chipDynoTransFactNoise.m: given a transcription factor, provides TFAs.
chipDynoActTransFactNoise.m: identifies significantly varying TFs.
chipStatMatrixInverterNoise.m: inverts block tridiagonal matrices for chipChip
chipDynoMaxDiff.m: computes most significant changes in TFAs
chipChipTextRead.m: reads TXT file for the Lee ChIP data files.
chipTuTextRead.m: assigns common names to probe IDs
demSpellmanp1p2.m: demonstrates dynamical chipCHIP on Spellman data.
chipDynoExpectationsFast.m: computes posterior expectations of TFA.
demSpellmanp5p2.m: demonstrates dynamical chipCHIP on Spellman data.
chipTextRead.m: reads TXT file for the Spellman data files.
chipDynoCompareX.m: compares connectivity for different pvals
chipStatMatrixInverter.m: inverts block tridiagonal matrices for chipChip
chipReduceVariables.m: reduce  number of variables in chipDyno model
chipDynoLikeStatGrad.m: gradient of CHIPDYNOLIKESTAT.
demSpellmanp2.m: demonstrates dynamical chipCHIP on Spellman data.
demSpellmanDynoStat.m: demonstrates dynamical chipCHIP on Spellman data.
chipDynoCompareTrans.m: compares list of active tfs at different p
chipDynoActTransFact.m: identifies significantly varying TFs.
chipDynoLoadDatap1.m: loads Spellman Data with Lee et al ChIP data.
demSpellmanp5.m: demonstrates dynamical chipCHIP on Spellman data.
chipDynoStatPostEstNoise.m: computes posterior expectations.
chipDynoLikeStat.m: marginal likelihood for chipChip dynamical model
chipDynoSignRels.m: finds the regulatory strength conserved rels
chipDynoLoadDatap1p2.m: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadDatap2p5.m: loads Spellman Data with Lee et al ChIP data.
chipDynoLoadData.m: loads Spellman Data with Lee et al ChIP data.
chipDynoGeneAct.m: given a gene, lists activators in decreasing order
demTuNoNoise.m: demonstrates dynamical chipCHIP on Tu data.
chipDynoLikeStatNoise.m: marginal likelihood for chipChip dynamical model
chipDynoGeneActNoise.m: given a gene, lists activators in decreasing order
chipDynoTransFact.m: provides gene-specific TFAs with errorbars.
demTu.m: demonstrates dynamical chipCHIP on Tu data.
chipDynoTuLoadData.m: loads Tu Data with Lee et al ChIP data.
chipDynoStatPostEst.m: computes posterior expectations 
chipDynoLoadDatap2.m: loads Spellman Data with Lee et al ChIP data.
chipDynoExpectationsFastNoise.m: computes posterior expectations TFA.
chipDynoNewTransFactNoise.m: tfs active for us and not for Tu et al.
chipDynoPerGenes.m: genes periodic according to Tu et al
chipDynoTransFactCorrel.m: correlations between transcription factors
chipDynoMultRegs.m: finds genes that have multiple regulators
chipDynoLikeStatNoiseGrad.m: gradient of CHIPDYNOLIKESTATNOISE.
chipDynoLoadDatap5.m: loads Spellman Data with Lee et al ChIP data.
