% Name        : [theVoter]=build_voter(nRows,nCols,centralHeight,theSlope)
% Description : Builds a voting matrix over the RR interval x sample number
%               space.
% Input       : nRows, nCols - Desired size of the voting matrix. Both
%               values must be odd. Otherwise, the voter might neither have
%               the exact desired size nor be symmetric.
%
%               centralHeight - Number of valued rows in the center. It
%               represents the allowed variation for the current
%               RR-interval.
%
%               theSlope - Slope for the allowed variation with time. A
%               slope of 1 means that for every time step (to the future or
%               to the past) the RR interval can change also one time step.
%               A value of .5 means that for every 2 time steps, the RR
%               interval can change 1 time step.
%
% Output      : theVoter - The voter matrix build as described next.
%
% Note        : The voter values not satisfying the constraints imposed by
%               centralHeight and theSlope are zero. Otherwise, their
%               values are proportional (within 0-1 range) to a Gaussian
%               PDF whose mean is at the center of the voter. The
%               covariance is set so that the voter covers the 4 sigma
%               interval.
% Author      : Antoni Burguera Burguera (2017)
%               antoni.burguera@uib.es
function [theVoter]=build_voter(nRows,nCols,centralHeight,theSlope)
    qRows=floor(nRows/2);
    qCols=floor(nCols/2);
    rStart=qRows-round(centralHeight/2)+1;
    q=zeros(qRows,qCols);
    for r=1:qRows
        for c=1:qCols
            if r>rStart-c*theSlope
                q(r,c)=mvnpdf([r;c],[qRows;1],[(nRows/8)^2 0;0 (nCols/8)^2]);
            end;
        end;
    end;
    q=q-min(min(q));
    q=q/max(max(q));
    vt=[q(:,end:-1:1) q(:,1) q];
    vb=[q(end:-1:1,end:-1:1) q(end:-1:1,1) q(end:-1:1,:)];
    vm=vt(end,:);
    theVoter=[vt;vm;vb];
return;