==============
=== CREDIT ===
==============

RRI correction in scenarios with intervals of extreme noise.

Author : Antoni Burguera
e-mail : antoni.burguera@uib.es

If you use this software, please cite the following paper:

Antoni Burguera.
Inter-Beat Interval Estimation from Extremely Noisy Single Lead Electrocardiograms
In Proceedings of the IFAC World Congress, Berlin - Virtual, 2020.

=========================
=== BUILD THE PACKAGE ===
=========================

Download and unpack the package, if necessary, wherever you want. Open Matlab and navigate to the package source directory. Type "ls" or "dir" and check that several .m files and one .c file are there. One of the .m files is "install.m".

The package is pure Matlab code except for the C file with MEX interface to Matlab. To build it, just execute the file "install.m":

>> install

or type

>> mex smooth_signal.c

Some warnings related to GCC version may appear, but no errors should occur. You should see "MEX completed successfully" appear once.

In case of errors, check your Matlab MEX configuration by typing "help mex" and following the instructions.

Also, smooth_signal is optional. In case you cannot compile it, just search where it is used and remove the function calls. The impact of smooth_signal on the final results is almost neglectable.

====================
=== RUN THE CODE ===
====================

Within the source directory there is a "main.m" file. Execute it to see an example:

>> main

This example loads a synthetic ECG generated with ECGSYN (https://www.physionet.org/physiotools/ecgsyn/Matlab/) corrupted with intervals of noise of different durations. By typing "dataSet.name" after the execution, some data about ECGSYN is displayed.

Even though "main" loads synthetic data, it shows the main usage of the software. In particular, the main function is:

[refinedRRI,refinedError,shortRRI,longRRI,initialRRI,guidanceRRI,guidanceError,voteTable]=process_rri(dataSet.Rdata(1,:),dataSet.Fs);

whose input parameters are the detected R peaks and the frequency, in Hz, at which the source ECG was sampled. The output includes both the refined RRI estimates and the original ones, as well as the beats classified as too short or too long.

Type

>> help process_rri

for more information. Also, it is advisable to navigate through the source code and read the comments and file headers to understand how the software works.

==========================
=== FULL USAGE EXAMPLE ===
==========================

To use the software with a Physionet ECG from Matlab, proceed as follows:

1) Install WFDB Matlab toolbox in case you don't have it (https://physionet.org/physiotools/matlab/wfdb-app-matlab/) and load one ECG

>> [tm,ecg,Fs]=rdsamp('nstdb/118e00',1);

In this example, the first parameter is the record path within WFDB search path. Type "help rdsamp" if it doesn't work or want to know which records are available.

The second parameter is the signal number. Some databases are multiparameter. Thus, make sure that the signal number corresponds to an ECG. A simple way to check it is to download all the record files and open the corresponding ".hea" file with a text editor.

2) Check Fs. In most cases, rdsamp does not place any value within Fs. If this happens, you must assign it manually. Again, check the ".hea" file (well documented in Physionet) to know the sampling frequency. Another option is to guess Fs from "tm", which is usually properly read, as follows:

>> Fs=round(length(tm)/(tm(end)-tm(1)));

3) Obtain the R-peak indexes using a QRS detector. Each detector has its own format. For example, you can use our QRS detector (https://github.com/aburguera/QRS_DETECT). If so, the R peaks are obtained as follows:

>> [Rdata,~]=process_ecg(ecg,Fs);

If you use this QRS detector, please refer to the README.TXT file in the corresponding github repository to learn how to cite us.

4) Execute the code

>> [refinedRRI,refinedError,shortRRI,longRRI,initialRRI,guidanceRRI,guidanceError,voteTable]=process_rri(Rdata(1,:),Fs);

5) That's it. Now you can plot and compare the refinedRRI and the initial RRI. Also, you can access the annotated beats for most databases within Physionet (see Physionet documentation). So, you can compare the initial RRI and the refined one with a ground truth.

Please, note that to properly plot/compare the data, you need a common time basis. See main.m to learn how to accomplish this.

Also, note that this description is Physionet and WFDB related. Actually, you just need a vector with the ECG samples (ecg), the sampling frequency (Fs) and a QRS detector.

=======================
=== TROUBLESHOOTING ===
=======================

This software package has been tested using Matlab 8.5.0.197613 (R2015a) running on Ubuntu 16.04LTS with GCC 5.4.0. Other configurations are untested.
