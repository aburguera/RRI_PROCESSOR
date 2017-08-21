% Name        : [initialRRI,guidanceRRI,guidanceError,voteTable]=build_guidance(Rind,Fs)
% Description : Builds the RRI guidance data from the detected R-peaks.
% Input       : Rind - Detected R-peaks. 1xM vector containing the ECG
%                       sample number where the R-peaks occured.
%
%               Fs - The frequency at which the original ECG was sampled
%                    (Hz)
%
% Output      : initialRRI - The RRI directly extracted from Rind.
%
%               guidanceRRI - The guidance RRI.
%
%               guidanceError - An error estimate for each guidance RRI
%
%               voteTable - The built voting table. See code.
%
% Author      : Antoni Burguera Burguera (2017)
%               antoni.burguera@uib.es
function [initialRRI,guidanceRRI,guidanceError,voteTable]=build_guidance(Rind,Fs)
    % ---
    % Pre-compute parameters, initialize storage and build vote table
    % ---

    % Clear space used for guidance nodes
    theNodes=[];
    % Compute the initial set of RR-intervals from the detected R-peaks.
    initialRRI=diff(Rind);
    % Initialize storage for the guidance error bounds.
    guidanceError=zeros(1,length(initialRRI));
    % The change in RR-interval allowed for the current beat.
    centralHeight=0.1*Fs;
    % The change in RR interval that is voted
    totalHeight=0.3*Fs;
    % Within how many beats this change can happen
    totalWidth=0.15*Fs;
    % Build the voter, ensuring dimensions are odd.
    theVoter=build_voter(round(totalHeight/2)*2+1,round(totalWidth/2)*2+1,centralHeight,2);
    % Build the vote table.
    voteTable=do_voting(initialRRI,theVoter);

    % ---
    % Build the guidance
    % ---

    % Get a threshold for voteTable by means of Otsu method.
    minVal=graythresh(voteTable);

    % Build the guidance node list. For each time step, search the most
    % voted RR interval. If the votes exceed the threshold, store it as a
    % node.
    for i=1:size(voteTable,2)
        [theValue,j]=max(voteTable(:,i));
        if theValue>minVal
            theNodes=[theNodes [i;j-1]];
        end;
    end;
    % Add a node at the beginning and the end if necessary.
    if theNodes(1,end)~=size(voteTable,2)
        theNodes=[theNodes [size(voteTable,2);theNodes(2,end)]];
    end;
    if theNodes(1,1)~=1
        theNodes=[[1;theNodes(2,1)] theNodes];
    end;

    % Build the guidance by interpolating and smoothing.
    guidanceRRI=interp1(theNodes(1,:),theNodes(2,:),1:length(initialRRI),'pchip');
    guidanceRRI=smooth_signal(guidanceRRI,round(totalWidth/2));

    % Let us assume that the obtained guidance (guidanceRRI) error is bound by
    % the perimeter of the binarized voting table.
    errorBounds=bwperim(im2bw(voteTable,minVal));

    lastError=30;   % Initialize to a reasonable value
    for i=1:size(errorBounds,2)
        % Search errorBounds for this time step
        j=find(errorBounds(:,i)~=0);
        % If they exist, set the error to the maximum difference between
        % the errorBounds and the current RR interval.
        if length(j)>0
            guidanceError(i)=max(abs(min(j)-guidanceRRI(i)),abs(max(j))-guidanceRRI(i));
            lastError=guidanceError(i);
        % If no errorBounds for this time step, use lastError and increase
        % it (uncertainty increases when no data available).
        else
            guidanceError(i)=lastError;
            lastError=lastError+.1;
        end;
    end;

    % Smooth the error estimate
    guidanceError=smooth_signal(guidanceError,round(totalWidth/2));

    % Complete the signals (diff outputs one item less)
    initialRRI=[initialRRI(1) initialRRI];
    guidanceRRI=[guidanceRRI(1) guidanceRRI];
    guidanceError=[guidanceError(1) guidanceError];
    voteTable=[voteTable(:,1) voteTable];
return;