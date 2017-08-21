% Load the ECG. This ECG is synthetic, generated with ECGSYN and corrupted
% with intervals of extreme noise. After loading, print the dataSet.name
% field for more information.
% The field dataSet.Rdata has been obtained from the QRS detector in:
% https://github.com/aburguera/QRS_DETECT
% by means of:
% [dataSet.Rdata,~]=process_ecg(dataSet.ecg,dataSet.Fs);
% Other QRS detector can be used. Refer to QRS_DETECT documentation for
% info about the output of process_ecg.
disp('* Loading data');
load syntheticDataCorrupted;

% Process the RRI.
disp('* Processing the RRI');
[refinedRRI,refinedError,shortRRI,longRRI,initialRRI,guidanceRRI,guidanceError,voteTable]=process_rri(dataSet.Rdata(1,:),dataSet.Fs);

% Plot the voting table and guidance.
disp('* Plotting');
figure;
image(voteTable*255);hold on;
plot(guidanceRRI,'k');hold on;
plot(guidanceRRI+guidanceError,'r');hold on;
plot(guidanceRRI-guidanceError,'r');
xlabel('Time (R-peak index)');
ylabel('RRI (samples)');
legend('Guidance','Guidance error interval','Location','SouthEast');
title('Voting table');
set(gca,'FontSize',14);
set(gca,'YDir','normal');

% Plot the output and ground truth.
theGT=[dataSet.annotation(1,2:end)/dataSet.Fs;diff(dataSet.annotation(1,:))/dataSet.Fs];
theRefined=[dataSet.Rdata(1,:)/dataSet.Fs;refinedRRI/dataSet.Fs];
theInitial=[dataSet.Rdata(1,:)/dataSet.Fs;initialRRI/dataSet.Fs];

figure;
plot(theInitial(1,:),theInitial(2,:),'k');hold on;
plot(theGT(1,:),theGT(2,:),'r');hold on;
plot(theRefined(1,:),theRefined(2,:),'g');hold on;
plot(theRefined(1,:),theRefined(2,:)+refinedError*(3/2)/dataSet.Fs,'b');hold on;
plot(theRefined(1,:),theRefined(2,:)-refinedError*(3/2)/dataSet.Fs,'b');
xlabel('Time (s)');
ylabel('RRI (s)');
xlim([theRefined(1,1) theRefined(1,end)]);
legend('Initial RRI','Gound truth','Processed RRI','3\sigma bound','Location','SouthEast');
title('RRI comparison');
set(gca,'FontSize',14);

% Plot too short and too long beats. As this example comes from a synthetic
% dataset, these outliers appear within the noisy areas. The non noisy
% areas almost do not have outliers.
figure;
plot((1:length(dataSet.ecg))/dataSet.Fs,dataSet.ecg,'k');hold on;
plot(dataSet.Rdata(1,shortRRI)/dataSet.Fs,dataSet.Rdata(2,shortRRI),'bo');hold on;
plot(dataSet.Rdata(1,longRRI)/dataSet.Fs,dataSet.Rdata(2,longRRI),'ro');
xlim([0 length(dataSet.ecg)/dataSet.Fs]);
xlabel('Time (s)');
ylabel('ECG amplitude');
legend('ECG','Detected short RRI','Detected long RRI');
title('Outlier detection');
set(gca,'FontSize',14);

disp('* Done');