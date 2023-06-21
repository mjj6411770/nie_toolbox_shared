
# Automated analysis of NIE signals 

It is a simple algorithm for using data-driven template matching in nano-impact electrochemistry signal analysis.





## Authors

- Ziwen Zhao [@ziwzh166](https://github.com/ziwzh166)
- Arunava Naha [@arun-naha](https://github.com/arun-naha)
- Sagar Ganguli [@gangulisagar](https://github.com/gangulisagar)
- Alina Sekretareva*[@alina-sekretareva](https://github.com/alina-sekretareva)


## Installation

The analysis algorithm is written in MATLAB. Please install the following add-ons to test the package: 

```matlab
    %% Signal processing toolbox
    %% Statistics and Machine Learning Toolbox
    %% Parallel Computing Toolbox
```
    
## The algorithm follows the flowchart below

![S1](https://github.com/ziwzh166/NIE_toolbox_shared/assets/100134089/32f6f62a-0d71-4637-b8f0-ac92df68d3da)




## To start in /test_simple
The directory "/test_simple" contains a streamlined version of the algorithm demonstrated on the NIE data for glucose electrooxidation on gold nanoparticles. The code is identical to the detailed one (in the "/test_detailed" directory), but it is more straightforward to use as many parameters are set by default. There are three .m files and two .txt files in this directory. "SignalAnalysis.m" is the algorithm that enables NIE data analysis,  "NIE_analysis_Simp.m" is the function-packed version of all the steps listed in the detailed version, "RUN_THIS_FILE.m" is the file that you can run to analyze your data with the developed algorithm. Additionally, there are two .txt files containing the NIE and blank data signals.
In the "RUN_THIS_FILE.m" 
there is only one line of the code:
```matlab
    [SpikeFeatures, SpikeLocation] = NIE_analysis_Simp("1a_AuNpsNIE.txt","1b_AuNpsBlank.txt","A","Oxi");
```
Inputs for the "RUN_THIS_FILE.m" in order from the left to the right are: NIE data file name, Blank data file name, the current unit in the txt file, and the reaction type, which can be "Oxi" or "Red" for oxidative or reductive reaction. The output is the features of the spikes and the location of the spikes in the original signal. The output is shown below:
The SpikeFeatures is a table with the following columns:
```matlab
    SpikeFeatures = 
    Slope_Left    Slope_Right    Time_Duration    Relative_Peak_Location    Height    Prominence    Area    Prominence_Duration_Ratio
```
![image](https://github.com/ziwzh166/NIE_toolbox_shared/assets/100134089/42aafab8-f71e-482d-926b-9bf301806780)

The SpikeLocation is a matrix with the following columns:
```matlab
    SpikeLocation = 
    Spike_Left   Spike_Right  Spike_Peak
```
![image](https://github.com/ziwzh166/NIE_toolbox_shared/assets/100134089/7af06c01-dab1-4ce7-9644-6debf9268798)

Note that the peak Location is expressed as Index instead of the time value in the original signal. The time value can be obtained by multiplying the index by the sampling rate. The sampling rate is 1/100 s in this example.

When you run the code, it will show two plots which are explained in "Generate the new templates marked spikes on an original signal related statistic" section below. The first plot is the original signal with the marked spikes, and the second plot is the statistic of the spikes. The statistic includes the number of spikes, the average height, the average area, and the average prominence of the spikes. 

## To start in /test_detailed

The directory "/test_detailed" contains two files in the .m format. The first one, "SignalAnalysis," is the algorithm that enables NIE data analysis, while the second .m file, "Display.mlx," is a live editor file that permits observing the results as they are generated. Additionally, there are four .txt files containing the NIE and blank signals discussed in the manuscript.

## Importing the data

```matlab
    Sig1 = SignalAnalysis(names_Str(1));
```

## Resampling 

This part of the code is only used for catalase NIE data and can be omitted for the data recorded with commonly used potentiostats.

```matlab
  Sig1 = Sig1.Preprocess("PA","N","Y","Y");
```

The first argument is used for changing the current unit. This simple toolbox converts all current scales into pA. The second argument is used to shift the time from 0, which can be input as 'yes' or 'no.' The third argument determines whether resampling is performed, and the last argument determines whether to show both the original and resampled data on a plot. The following plot shows the result of resampling:
![Resampling](https://user-images.githubusercontent.com/100134089/224944121-91084fef-a1f0-4e92-b4d2-900f48043e30.svg)

## Signal denoising 

As stated in the paper, we denoise the signal by using a lowpass filter and understanding the noise frequency by the stft function. The denoise method has two inputs. The first one is the cutoff frequency, and the second is the order of the filter to apply. In general, the higher the order, the stronger the filter, but the longer the signal delay. If you call the method of the class by leaving the first argument blank, as shown below:
```matlab
  Sig1 = Sig1.Denoise([],10);
```
The output will be the 3D plot below. The time axis may differ:

![stftCATA1](https://user-images.githubusercontent.com/100134089/234817164-388d9127-ccb7-436e-9be0-7558cd55fafa.png)


If you set the cutoff frequency:

```matlab
  Sig1 = Sig1.Denoise(15,10);
```
It will output the following curve:
![DenoiseCATA](https://user-images.githubusercontent.com/100134089/224957690-a83e87b4-bb0b-4373-a720-10b8a3be4b18.svg)

It will leave an input window on the command window asking the user if they want to change the order of the filter and the possible order values.

## Background trend removal 
The background subtraction is made by making the smooth data using the rloess. The method has two inputs: the window size and the 'Y' or 'N' argument to select whether to show the plot for the fitted line and offset signal. In general, a larger window size results in a better-fitting trend line.    

```matlab
  Winsize = 100;
  Sig1 = Sig1.BgSub(Winsize,'y');
  Back1 = Back1.BgSub(1000,'N');
```
![BgSubCat1](https://user-images.githubusercontent.com/100134089/224961351-3191704a-ea74-41f8-b866-184037d515aa.svg)

## Initial spike sampling via the conventional heigh threshold method 
This step uses the blank signal to determine the height threshold. The sampling is then performed by this code:

```matlab
Sig1 = Sig1.FlipFindPeak(Back1.Test_signal_offset,'Red','Y',"Offset"); 
```
This method accepts four arguments (in the order from the left to the right): the de-trended blank signal, 'Red' or 'Oxi' showing whether the NIE signal is oxidative or reductive, 'Y' or 'N' to choose whether to show the found peaks, 'Offset' or 'Original' to choose whether to show the found peaks on the offset or original trend signal. The example of the output:
![FindPeaksCAT](https://user-images.githubusercontent.com/100134089/224968188-8bab673a-283b-4c9e-990f-4a18a26fb868.svg)

## Initial spike feature extraction and automated spike grouping 

To identify different spike-shape templates the K-means clustering is performed on the sampled spikes. The method is as 

```matlab
Sig1 = Sig1.GeRawTrainSet('Y','Y');
```
The first argument 'Y' or 'N' is to choose whether to use or not the elbow method, and the second 'Y' or 'N' is to choose whether to use or not the silhouette score. 
#### Elbow method:
![ElbowCAT](https://user-images.githubusercontent.com/100134089/224988058-98917542-1926-47d3-aefd-f1bb2a75ce32.svg)


#### Silhouette Score:
![SC_CAT](https://user-images.githubusercontent.com/100134089/224986902-edcbbe18-0032-42e8-a6cf-07d80d376b41.svg)


## Template generation
The next step is to generate the representative templates based on the K-means results. The spikes belonging to the same cluster are averaged to obtain the raw templates using the following command:
```matlab
Sig1 = Sig1.KmeansGeRawSigTem(3,'Y',"Y");
```
it will show the raw templates and a bar chart showing % composition of the signal by different templates
![TemRawCAT](https://user-images.githubusercontent.com/100134089/224986136-49ab4ef3-6b7f-46ba-88df-ae843e7fe360.svg)
![BarCAT](https://user-images.githubusercontent.com/100134089/224989387-41b8a0a0-6ba9-401f-b0e7-9ffc79e0acf9.svg)

## Template tuning
There are some templates which are templates for the noisy spikes, and some templates need to regulate two sides to avoid involving too much background trend; the first black one is removed by the method output as the last argument works as whether to plot 
```matlab
Sig1 = Sig1.RawtemplatesReguFunc([2,3],'Y');
```
![Regutem](https://user-images.githubusercontent.com/100134089/224990028-70f0515b-8e6c-48ab-9137-b300eacfd2d7.svg)

## Template Matching 
templates matching is performed by NCC coefficient the details can be found in the paper,the NCC coefficient shows a cosine similarity, given by the method:
```matlab
Sig1 = Sig1.Templatematching('Y');
```
the argument decides whether to show the plot as below, xy axis is linked 
![SimCurve](https://user-images.githubusercontent.com/100134089/224991829-1b20b68b-26f0-4990-b84c-8269cdb2fa37.svg)

## Numerical filtering and interval merging 
The matching is given the matched spikes with the noisy scale, which can be filtered by the numerical filter given by the method below:
```matlab
Sim = 0.9;
StdCoeff = 0.35;
HeightWidth = 0.35;
Sig1 = Sig1.TemplatematchingFiltering(Sim,StdCoeff,HeightWidth,'Y','Y');
```
The sim is the similarity filter; the StdCoeff and HeightWidth compare each matched spike with the template corresponding value times coefficients. Followed by two plot arguments; first one show the matched interval by two different templates, and the second shows the merged different matched interval by different colors:
![MatchedTeminte](https://user-images.githubusercontent.com/100134089/224993990-5aa9a6f1-ed88-43d6-808c-d4ea28757ad9.svg)
![MatchedMerInte](https://user-images.githubusercontent.com/100134089/224994486-4957a062-d34b-45d3-92b2-a4b0d403ec5e.svg)

## Physical info extracted and regenerated the templates
With two well-defined sides, we can extract the information as what we want to have. There are eight parameters get extracted. They can be found in the paper
these are given by the code: 
```matlab
Sig1 = Sig1.GeAMTrainSet('Y','Y',"N");
```
The first argument is in case the numerical filtering is still not good enough, it will take the height threshold that calculates from the blank to filter further the peak, and the following two are the plot relates to select numbers to K-means we only show the elbow methods here for the other one you can run the test by your own.
![Elbow2](https://user-images.githubusercontent.com/100134089/225002860-d110a05c-4e00-4e59-a475-ecd4df06687f.svg)

## Generate the new templates marked spikes on an original signal related statistic
after reclustering the physical information we can regenerate the templates as the representatives for the signal spikes by the given method 
```matlab
Sig1 = Sig1.KmeansGeAMSigTem(2,"Y","Y");
```
the first one is the centroid numbers; the second shows the found spikes on the original denoised signal and some related statistic 
#### Marked original signal:
![CAT_Marked](https://user-images.githubusercontent.com/100134089/234838710-e7a32470-9484-48bf-aee5-08708cf6b71c.png)
#### Reclustering and showing the statistic:
![F5](https://user-images.githubusercontent.com/100134089/234835731-3987d5ae-27ee-43fe-8dc4-9b2836441469.png)
