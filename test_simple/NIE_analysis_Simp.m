function [Features,SpikeIndex] = NIE_analysis_Simp(Datafilename,Blankfile,unit,ReactionStates)
%[Features,SpikeIndex] = NIE_analysis_Simp(Datafilename,unit)
%   just a function to compress all the steps 
%   Datafilename: the name of the data file
%   Blankfile: the name of the blank file
%   unit: the unit of the current by default is "A"
%   ReactionStates: the reaction states of the data file as "Oxi" or "Red" by default is "Oxi"
%   Features: the features of the spikes
%   SpikeIndex: the index of the spikes contains left side, right side and peak index
%   for details check the folder test_detailed
if nargin < 3
    unit = "A";
    ReactionStates = "Oxi";
end
if nargin < 4
    ReactionStates = "Oxi";
end

% Input the data file and blank file
Sig = SignalAnalysis(Datafilename);
Blank = SignalAnalysis(Blankfile);
% Preprocess the data file and blank file
Sig = Sig.Preprocess(unit,"Y","N","N");
Blank = Blank.Preprocess(unit,"Y","N","N");
% Denoise the data file
Sig = Sig.Denoise(15,10,"No");
% Baseline correction
Sig = Sig.BgSub;
Blank = Blank.BgSub;
% First spike sampling
Sig = Sig.FlipFindPeak(Blank.Test_signal_offset,ReactionStates,"No");
% First feature extraction
Sig = Sig.GeRawTrainSet('No','No');
% First clustering
Sig = Sig.KmeansGeRawSigTem;
% Template tuning
Sig = Sig.RawtemplatesReguFunc([],'N');
% Template matching 
Sig = Sig.Templatematching('No');
% Numerical filtering for second spike sampling
Sig = Sig.TemplatematchingFiltering;
% Second feature extraction
Sig = Sig.GeAMTrainSet('N','N',"N");
% Second clustering
Sig = Sig.KmeansGeAMSigTem([],"Y","Y");
% Extracted features and spike index
Features = Sig.AMtrainingSet;
SpikeIndex = Sig.TemplatateMatchedInte;






end

