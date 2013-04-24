function [geneName, annotation, data] = chipChipTextRead(file1, file2)

% CHIPCHIPTEXTREAD reads TXT file for the Lee ChIP data files.
% CHIPDYNO toolbox
% chipChipTextRead.m version 1.5
% FORMAT function [geneName, annotation, data] = chipChipTextRead(file1, file2)
% DESC reads TXT file for the Lee ChIP data files.
% ARG file1 : the geneName and annotation file
% ARG file2 : the data file
% RETURN geneName: geneNames
% RETURN annotation: annotation of the geneNames
% RETURN data : TFA at different experimental point
% COPYRIGHT : Neil D. Lawrence, 2005
% COPYRIGHT : Guido Sanguinetti, 2005
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipTextRead, chipTuTextRead

[geneName,annotation] = ...
    textread(file1,'%q %q %*[^\n]',...
             'headerlines',2,'whitespace','','delimiter','\t');

data = load(file2);
