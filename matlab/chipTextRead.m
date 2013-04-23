function [geneName, data] = chipTextRead(file)

% CHIPTEXTREAD reads TXT file for the Spellman data files.
% CHIPDYNO toolbox
% chipTextRead.m version 1.5
% FORMAT [geneName, data] = chipTextRead(file)
% DESC reads TXT file for the Spellman data files.
% ARG file : the Spellman data file
% RETURN geneNames : Gene names
% RETURN data : point estimate of the expression level
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipTuTextRead


% file ia a string containing the file name and the extension.
% data is a matrix with 24 columns and N( number of genes) rows   
[geneName,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14...
    x15,x16,x17,x18,x19,x20,x21,x22,x23,x24]=...
    textread(file,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'headerlines',1,'whitespace','','delimiter','\t');

data=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24];
