
# NIE data-driven template matching 

It's a simple toolbox for using the data-driven template matching in nano-impact electrochemistry signal analysis




## Authors

- Ziwen Zhao [@ziwzh166](https://github.com/ziwzh166)
- Arunava Naha [@arun-naha](https://github.com/arun-naha)
- Sagar Ganguli
- Alina Sekretareva*


## Installation

the project is based on Matlab, please install the following add-ons to test the package 

```matlab
    %% Signal processing toolbox
    %% Statistics and Machine Learning Toolbox
    %% Parallel Computing Toolbox
```
    
## This toolbox will follow the flow chart below

![Flow chart](https://user-images.githubusercontent.com/100134089/224944354-9ec909b2-6663-45c2-b86c-1dfd22031aa9.svg)



###
## To start

In the folder /test, you can find two .m files, SignalAnalysis is the toolbox for performing the analysis, and Display.m is the live editor file in which you can see the results step by step. Furthermore, the other four .txt files are the NIE and blank signals we mentioned in the paper.

In the first section are created to import the data and create new classes
which is given by the first name string in the folder and is imported into the class.

```matlab
    Sig1 = SignalAnalysis(names_Str(1));
```

## Resampling 

As we state in the paper, the unevenly sampled signal is integrated into the preprocess method.

```matlab
  Sig1 = Sig1.Preprocess("PA","N","Y","Y");
```

The first argument is made for changing the current unit. This simple toolbox will convert all scales of current into pA; the second argument is whether or not to shift the time from 0, which can be input as 'yes' or 'no.' Then the third argument is whether the resample is performed or not, and the last one is if to show the plot with an original and resampled signal. 
It will be the plot as below: 
![Resampling](https://user-images.githubusercontent.com/100134089/224944121-91084fef-a1f0-4e92-b4d2-900f48043e30.svg)

## Denoise 

As we state in the paper, we denoise the signal by the lowpass filter and understand the noise frequency by stft function, in the method we integrated,
there are two inputs, fisrt one is the frequecy to cut off and second is the orders of the filter to apple, in general the higher the order the stronger the filetr it will be but stonger delay to the signal. If you call the method of class by leaving first argument blank as below 

```matlab
  Sig1 = Sig1.Denoise([],10);
```
It will output the 3D plot as below, the time axis may differ:

![stftCATA1](https://user-images.githubusercontent.com/100134089/224954832-f8a181a1-0020-407c-a0c7-f6e5dad1f0d2.svg)

In case you set the frequency like below: 

```matlab
  Sig1 = Sig1.Denoise(15,10);
```
it will ouput the curve as one below:
![DenoiseCATA](https://user-images.githubusercontent.com/100134089/224957690-a83e87b4-bb0b-4373-a720-10b8a3be4b18.svg)

it will leave a input window on the command window to ask if the user would like to change the order and the orders to be changed.

## Background Subtraction 
The background suntraction is made by making the smoothdata based on the roloess, you can input the argument as, the bigger the window size the better fitting the trend line, the last argument is wehter or not to show the plot for the fitted line and offset signal, the 

```matlab
  Winsize = 100;
  Sig1 = Sig1.BgSub(Winsize,'y');
  Back1 = Back1.BgSub(1000,'N');
```
![BgSubCat1](https://user-images.githubusercontent.com/100134089/224961351-3191704a-ea74-41f8-b866-184037d515aa.svg)

## Sampling some peaks based on height threshold 
This step is based on the blank signal to generate the height threshold and perform sampling follow by this code:

```matlab
Sig1 = Sig1.FlipFindPeak(Back1.Test_signal_offset,'Red','Y',"Offset"); 
```
this method accept four arguments, the first one is the detrended blank signal, the next one is either the NIE signal is oxidative or reductive. in 'Red' or 'Oxi', the the third argument is if to show the plot and show the found peaks on the offset or original tren signal, the output as below:
![FindPeaksCAT](https://user-images.githubusercontent.com/100134089/224968188-8bab673a-283b-4c9e-990f-4a18a26fb868.svg)

## Generating the traning sets and Clustering 

Some parameters are extacted to distinguish the types of the spikes and noise spikes, performed by the K-means, the method is as 

```matlab
Sig1 = Sig1.GeRawTrainSet('Y','Y');
```
the first argument is based on elbow method and the second is by silhouette score, input 'y' or 'n' for the corresponding methods 
#### Elbow method:
![ElbowCAT](https://user-images.githubusercontent.com/100134089/224988058-98917542-1926-47d3-aefd-f1bb2a75ce32.svg)


#### Silhouette Score:
![SC_CAT](https://user-images.githubusercontent.com/100134089/224986902-edcbbe18-0032-42e8-a6cf-07d80d376b41.svg)


## Templates generation
The next step is to generate the representative templates based on the K-means results, the spikes belonging to same cluster are averaged, obtain the raw templates; the method is made by the following command:
```matlab
Sig1 = Sig1.KmeansGeRawSigTem(3,'Y',"Y");
```
it will show the raw templates and bar chart for the sum of different template
![TemRawCAT](https://user-images.githubusercontent.com/100134089/224986136-49ab4ef3-6b7f-46ba-88df-ae843e7fe360.svg)
![BarCAT](https://user-images.githubusercontent.com/100134089/224989387-41b8a0a0-6ba9-401f-b0e7-9ffc79e0acf9.svg)

## Templates Regulation
There are some templates which are the templates for the noisy spikes, and some templates need to regulate two sides to avoid involving two much background trend, the first black one is removed by the method output as, the last argument works as whether to plot 
```matlab
Sig1 = Sig1.RawtemplatesReguFunc([2,3],'Y');
```
![Regutem](https://user-images.githubusercontent.com/100134089/224990028-70f0515b-8e6c-48ab-9137-b300eacfd2d7.svg)

## Templates Matching 
templates matching is performed by NCC coefficient the deatils can be found in paper,the NCC coefficient shows a cosine similarity, given by the method:
```matlab
Sig1 = Sig1.Templatematching('Y');
```
the argument decides whetherto show the plot as below, xy axis are linked 
![SimCurve](https://user-images.githubusercontent.com/100134089/224991829-1b20b68b-26f0-4990-b84c-8269cdb2fa37.svg)

## Numerical filtering and interval mergering 
The matching given the matched spikes with the noisy scale, which can be filtered by the numercial filter, given by the method below:
```matlab
Sim = 0.9;
StdCoeff = 0.35;
HeightWidth = 0.35;
Sig1 = Sig1.TemplatematchingFiltering(Sim,StdCoeff,HeightWidth,'Y','Y');
```
the sim is the similarity filter, the StdCoeff and HeightWidth compare each matched spike with the templates corresponding value times coefficients. Following by two plot arguments, first one show the matched interval by two different templates, and the second shows the merged different matched interval by different colors:
![MatchedTeminte](https://user-images.githubusercontent.com/100134089/224993990-5aa9a6f1-ed88-43d6-808c-d4ea28757ad9.svg)
![MatchedMerInte](https://user-images.githubusercontent.com/100134089/224994486-4957a062-d34b-45d3-92b2-a4b0d403ec5e.svg)

## Physcial info extracted and regenerating the templates
With well defined two sides we can extarct the information as what we want to have there are eigth parameters get extracted. They can be found in paper
these are given by the code: 
```matlab
Sig1 = Sig1.GeAMTrainSet('Y','Y',"N");
```
the first arguemtns is in case, the numerical filetring is still not good enough, it will take the height threshold that calculate from the blank to filter further the peak and next two are the plot relates to select numbers to K-means we only show the elbow methods here for the other one you can run the trest by yourself
![Elbow2](https://user-images.githubusercontent.com/100134089/225002860-d110a05c-4e00-4e59-a475-ecd4df06687f.svg)

## Generate the new templates see it's on the original signal and some statistic data
![Recluster](https://user-images.githubusercontent.com/100134089/225004286-c07e9bc2-ebbb-4506-85ab-c75cb30980ce.svg)
![STA](https://user-images.githubusercontent.com/100134089/225004251-b02828cc-76f6-4936-a4ae-5bfbdaacf8c5.svg)
