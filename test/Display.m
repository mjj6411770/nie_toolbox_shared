%% ls1
dinfo = dir('*.txt');
names_cell = {dinfo.name}';
names_Str = convertCharsToStrings(names_cell)
%% import data
Sig1 = SignalAnalysis(names_Str(2));
Sig2 = SignalAnalysis(names_Str(3));
Back1 = SignalAnalysis(names_Str(1));
Back2 = SignalAnalysis(names_Str(4));
%% preprocess
Sig1 = Sig1.Preprocess("PA","N","Y","Y");
Back1 = Back1.Preprocess("PA","N","Y","Y");
Sig2 = Sig2.Preprocess("A","Y","N","N");
Back2 = Back2.Preprocess("A","Y","N","N");
%% Denoise
Sig1 = Sig1.Denoise(15,10);
%%
Sig2 = Sig2.Denoise(10,10);
%% BackgroundSubtraction Sig1
Winsize = 100;
Sig1 = Sig1.BgSub(Winsize,'y');
Back1 = Back1.BgSub(1000,'N');
%% BackGroundsub Sig2
Winsize = 600;
Sig2 = Sig2.BgSub(Winsize,'y');
Back2 = Back2.BgSub(1000,'N');
%% flip find peaks
Sig1 = Sig1.FlipFindPeak(Back1.Test_signal_offset,'Red','Y',"Offset");
Sig2 = Sig2.FlipFindPeak(Back2.Test_signal_offset,'Oxi','Y',"Offset");
%% TraningSet generation
Sig1 = Sig1.GeRawTrainSet('Y','Y');
Sig2 = Sig2.GeRawTrainSet('Y','Y');
%% Clustering and Raw unregularated templates generation
Sig1 = Sig1.KmeansGeRawSigTem(3,'Y',"Y");
Sig2 = Sig2.KmeansGeRawSigTem(5,'Y',"Y");
%% Regulated Templates
Sig1 = Sig1.RawtemplatesReguFunc([1,2],'Y');
Sig2 = Sig2.RawtemplatesReguFunc([1,2,4,5],'Y');
%% Template Matching
Sig1 = Sig1.Templatematching('Y');
Sig2 = Sig2.Templatematching('Y');
%% Template filtered the peaks Sig1
Sim = 0.9;
StdCoeff = 0.35;
HeightWidth = 0.35;
Sig1 = Sig1.TemplatematchingFiltering(Sim,StdCoeff,HeightWidth,'Y','Y');
%% Template filtered the peaks Sig2
Sim = 0.9;
StdCoeff = 0.35;
HeightWidth = 0.35;
Sig2 = Sig2.TemplatematchingFiltering(Sim,StdCoeff,HeightWidth,'Y','Y');
%% TraninignSet ReClustering
Sig1 = Sig1.GeAMTrainSet('Y','Y',"N");
Sig2 = Sig2.GeAMTrainSet('N','Y',"N");
%% ReClustering and generate the templates
Sig1 = Sig1.KmeansGeAMSigTem(2,"Y","Y");
Sig2 = Sig2.KmeansGeAMSigTem(3,"Y","Y");
%%
Sig3 = Sig2.KmeansGeAMSigTem(3,"Y","Y");
%% Rematching Once 
% it will replace the regulated raw templates by the new templates 
% Sig1 = Sig1.Rematching(6,"Y","Y");
Sig2 = Sig2.Rematching(5,"Y","Y");
%% Compare with the conventional method
Sig1 = Sig1.ConGaussian('Y');
Sig2 = Sig2.ConGaussian('Y');
%% for rematching
Sig3 = Sig1;
Sig4 = Sig2;
%% rematching 
Sig3 = Sig3.Rematching(5,"Y");
Sig4 = Sig4.Rematching(5,"Y");