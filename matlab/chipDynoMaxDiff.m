function f=chipDynoMaxDiff(TF,TFErrorDiff);

% CHIPDYNOMAXDIFF computes most significant changes in TFAs
% CHIPDYNO toolbox
% chipDynoMaxDiff.m version 1.2
% FORMAT f=chipDynoMaxDiff(TF,TFErrorDiff);
% DESC computes most significant changes in TFAs
% ARG TF: gene specific transcription factor activity
% ARG TFErrorDiff: error in gene specific transcription factor activity
% RETURN f : most significant changes in TFAs
% COPYRIGHT : Neil D. Lawrence, 2006
% COPYRIGHT : Guido Sanguinetti, 2006
% MODIFICATIONS : Muhammad A. Rahman, 2013
% SEEALSO : chipDynoExpectationsFast, chipDynoExpectationsFastNoise

nTargets=size(TF,1);
npts=size(TF,2);
preDiffs=zeros(npts,npts);
diffs=zeros(npts,npts,nTargets);
f=zeros(1,nTargets);
for i=1:nTargets
    for j=2:npts-1
        for l=j:npts
            preDiffs(j,l)=TF(i,j)-TF(i,l);
        end
    end
    diffs(:,:,i)=preDiffs-preDiffs';
    f(i)=max(max(diffs(:,:,i)./TFErrorDiff(:,:,i)));
end
