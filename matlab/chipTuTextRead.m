function [ORF, data, vars]=chipTuTextRead()

% CHIPTUTEXTREAD assigns common names to probe IDs
% CHIPDYNO toolbox
% chipTuTextRead.m version 1.5
% FORMAT [ORF, data, vars]=chipTuTextRead()
% DESC assigns common names to probe IDs
% RETURN ORF : Gene names
% RETURN data : point estimate of the expression level
% RETURN vars : uncertainty of the expression level
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipTextRead


[IDs,crap,ORF,geneCN]=textread('./data/MetabolData/dictionary.txt','%q%q%q%q');
IDSnew=textread('./data/MetabolData/probeIDTu.txt','%q');
data=load('./data/MetabolData/YeastMetabolism_exprs.txt');
vars=load('./data/MetabolData/YeastMetabolism_se.txt');
