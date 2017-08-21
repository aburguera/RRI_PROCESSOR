% Name        : [voteTable]=do_voting(theBeats,theVoter)
% Description : Builds the vote table in the space RR-interval x sample
%               number.
% Input       : theBeats - 1xN vector of RR-intervals.
%
%               theVoter - Voting matrix. Ideally constructed with
%               build_voter. The number of rows and columns MUST be ODD.
%
% Output      : voteTable - The resulting vote table, built as described
%               next.
%
% Note        : The rows in the vote table correspond to RR-interval
%               values, ranging from 1 to the maximum number of beats
%               within theBeats. The columns correspond to the time step,
%               ranging from 1 to the number of theBeats.
%               The voter's center is placed at each beat in theBeats and
%               its values are summed in the votes table. At the end, the
%               values are emphasized (squared) and normalized (between 0
%               and 1).
% Note        : The size of the vote table could be improved to:
%               - Avoid losing time when a few wrong extra-large beats lead
%               to large number of rows.
%               - Allow beats larger than the ones imposed by the maximum.
% Author      : Antoni Burguera Burguera (2017)
%               antoni.burguera@uib.es
function [voteTable]=do_voting(theBeats,theVoter)
    % Pre-compute parameters and initialize storage
    nBeats=length(theBeats);
    voteTable=zeros(max(theBeats),nBeats);
    if ((mod(size(theVoter,1),2)==0) || (mod(size(theVoter,2),2)==0))
        error('[ERROR] theVoter must have an odd number of rows and columns');
    end;
    halfHeight=(size(theVoter,1)-1)/2;
    halfWidth=(size(theVoter,2)-1)/2;
    % Build the vote table
    for i=1:nBeats
        % Compute boundaries of voter and vote table for the cases in which
        % the voter does not fully lie inside the vote table.
        bl=max(1,i-halfWidth);
        br=min(nBeats,i+halfWidth);
        bu=min(max(theBeats),theBeats(i)+halfHeight);
        bd=max(1,theBeats(i)-halfHeight);
        svl=bl-(i-halfWidth)+1;
        svd=bd-(theBeats(i)-halfHeight)+1;
        svr=2*halfWidth-((i+halfWidth)-br)+1;
        svu=2*halfHeight-((theBeats(i)+halfHeight)-bu)+1;
        % Update the vote table
        voteTable(bd:bu,bl:br)=voteTable(bd:bu,bl:br)+theVoter(svd:svu,svl:svr);
    end;
    % Square the vote table to emphasize the differences.
    voteTable=voteTable.^2;
    % Normalize it to values between 0 and 1
    voteTable=voteTable-min(min(voteTable));
    voteTable=voteTable/max(max(voteTable));
return;

