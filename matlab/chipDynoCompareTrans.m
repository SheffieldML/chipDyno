function [X,Y]=chipDynoCompareTrans(actTrans1,actTrans2)

% CHIPDYNOCOMPARETRANS compares list of active tfs at different p
%
%	Description:
%	[X,Y]=chipDynoCompareTrans(actTrans1,actTrans2)
%% 	chipDynoCompareTrans.m version 1.1

X=[];
for i=1:size(actTrans1,1)
  X=[X,sum((strcmp(actTrans1(i),actTrans2)))];
end
Y=[];
for i=1:size(actTrans2,1)
  Y=[Y,sum((strcmp(actTrans2(i),actTrans1)))];
end
