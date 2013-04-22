function [geneName, data] = chipTextRead(file)

% CHIPTEXTREAD reads TXT file for the Spellman data files.
%
%	Description:
%	[geneName, data] = chipTextRead(file)
%% 	chipTextRead.m version 1.5


% file ia a string containing the file name and the extension.
% data is a matrix with 24 columns and N( number of genes) rows   
[geneName,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14...
    x15,x16,x17,x18,x19,x20,x21,x22,x23,x24]=...
    textread(file,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'headerlines',1,'whitespace','','delimiter','\t');

data=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24];
