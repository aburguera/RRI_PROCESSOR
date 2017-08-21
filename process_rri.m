% Name        : [refinedRRI,refinedError,shortRRI,longRRI,initialRRI,guidanceRRI,guidanceError,voteTable]=process_rri(Rind,Fs)
% Description : Processes the RRI information in order to improve it.
% Input       : Rind - Detected R-peaks. 1xM vector containing the ECG
%                       sample number where the R-peaks occured.
%
%               Fs - The frequency at which the original ECG was sampled
%                    (Hz)
%
% Output      : refinedRRI - The processed RRI, ideally improving the
%               initial one. It is a 1xN vector of RRI measured in time
%               steps. In particular, refinedRRI(i) informs about the
%               length from the RRI that ended at time Rind(i). Note
%               that it may NOT coincide with Rind(i)-Rind(i-1), as
%               this is a refined version.
%
%               refinedError - The 2Sigma value for each item in
%               refinedRRI.
%
%               shortRRI - Indexes relative to Rind of the beats detected
%               as too short.
%
%               longRRI - Indexes relative to Rind of the beats detected
%               as too long.
%
%               initialRRI - The RRI directly extracted from Rind. See
%               build_guidance.
%
%               guidanceRRI - The guidance RRI. See build_guidance.
%
%               guidanceError - An error estimate for each guidance RRI.
%               See build_guidance.
%
%               voteTable - The built voting table. See build_guidance.
%
% Author      : Antoni Burguera Burguera (2017)
%               antoni.burguera@uib.es
function [refinedRRI,refinedError,shortRRI,longRRI,initialRRI,guidanceRRI,guidanceError,voteTable]=process_rri(Rind,Fs)
    % Build the guidance and get the initial estimates.
    [initialRRI,guidanceRRI,guidanceError,voteTable]=build_guidance(Rind,Fs);

    % Pre-compute some parameters
    Q=((0.2*Fs)/2)^2;  % Prediction model noise. Interpretation: 95% of cases
                        % (i.e. 2 sigma) RRI will be 0.2 s larger or shorter
                        % than the previous ones.

    % Initialize state vector and variance
    X=initialRRI(1);
    P=0;

    % Initialize storage
    refinedRRI=X;
    refinedError=P;
    shortRRI=[];
    longRRI=[];

    % Loop through all the data
    for i=2:length(initialRRI)
        % ---
        % KF prediction
        % ---

        % During prediction, the last RRI estimate (X) is not changed. Only the
        % variance is updated to take uncertainty into account.
        P=P+Q;

        % ---
        % KF update
        % ---

        % Store short and long beats by comparing them with the 2sigma bound.
        if initialRRI(i)<(X-2*sqrt(P))
            shortRRI=[shortRRI i];
        elseif initialRRI(i)>(X+2*sqrt(P))
            longRRI=[longRRI i];
        end;

        % Build the measurement vector, including the initial estimates and the
        % guidance.
        Z=[initialRRI(i);guidanceRRI(i)];
        H=[1;1];
        R=[Q 0; 0 (guidanceError(i)/2)^2];

        % Execute the KF update
        y=Z-H*X;
        S=H*P*H'+R;
        K=P*H'*inv(S);
        X=X+K*y;
        P=(eye(size(P))-K*H)*P;

        % Store the results
        refinedRRI=[refinedRRI X];
        refinedError=[refinedError 2*sqrt(P)];
    end;
return;