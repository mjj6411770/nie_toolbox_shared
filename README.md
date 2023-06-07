
# Automated analysis of NIE signals 

It is a simple algorithm for using data-driven template matching in nano-impact electrochemistry signal analysis.1




## Authors

- Ziwen Zhao [@ziwzh166](https://github.com/ziwzh166)
- Arunava Naha [@arun-naha](https://github.com/arun-naha)
- Sagar Ganguli [@gangulisagar](https://github.com/gangulisagar)
- Alina Sekretareva*[@alina-sekretareva](https://github.com/alina-sekretareva)


## Installation

The analysis algorithm is based on Matlab. Please install the following add-ons to test the package: 

```matlab
    %% Signal processing toolbox
    %% Statistics and Machine Learning Toolbox
    %% Parallel Computing Toolbox
```
    
## The algorithm follows the flowchart below

![S1b](https://github.com/ziwzh166/NIE_toolbox_shared/assets/100134089/9b5d852f-9b35-4faf-9247-8da18381d116)
com/ziwzh166/NIE_toolbox_shared/assets/100134089/4b97ec58-aadb-4ac4-a5e3-ef5c66c4b943)




###
## To start

Inside the "/test_detailed" directory, there are two files in the .m format. The first one, "SignalAnalysis," is the algorithm that enables NIE data analysis, while the second .m file, "Display.mlx," is a live editor file that permits to observe the results as they are generated. Additionally, there are four .txt files containing the NIE and blank signals discussed in the manuscript.
Inside the "/test_simple" it's the packed example with the glucose photo oxidation example. The code is the same as the detailed one, but it's more straightforward to use.
## Importing the data

```matlab
    Sig1 = SignalAnalysis(names_Str(1));
```

## Resampling 

This part of the code is only used for catalase NIE data and can be omited for the data recorded with commonly used potentiostats.

```matlab
  Sig1 = Sig1.Preprocess("PA","N","Y","Y");
```

The first argument is used for changing the current unit. This simple toolbox converts all current scales into pA. The second argument is used to shift the time from 0, which can be input as 'yes' or 'no.' The third argument determines whether resampling is performed, and the last argument determines whether to show the plot with an original and resampled signal. The following plot shows the result of resampling:
![Resampling](https://user-images.githubusercontent.com/100134089/224944121-91084fef-a1f0-4e92-b4d2-900f48043e30.svg)

## Denoising 

As stated in the paper, we denoise the signal by using a lowpass filter and understanding the noise frequency by the stft function. The method we integrated has two inputs. The first one is the frequency to cut off, and the second is the order of the filter to apply. In general, the higher the order, the stronger the filter, but the stronger the delay to the signal. If you call the method of the class by leaving the first argument blank, as shown below:
```matlab
  Sig1 = Sig1.Denoise([],10);
```
The output will be the 3D plot below. The time axis may differ:

![stftCATA1](https://user-images.githubusercontent.com/100134089/234817164-388d9127-ccb7-436e-9be0-7558cd55fafa.png)


If you set the frequency like this:

```matlab
  Sig1 = Sig1.Denoise(15,10);
```
It will output the following curve:
![DenoiseCATA](https://user-images.githubusercontent.com/100134089/224957690-a83e87b4-bb0b-4373-a720-10b8a3be4b18.svg)

It will leave an input window on the command window to ask if the user would like to change the order and the orders to be changed.

## Background Subtraction 
The background subtraction is made by making the smooth data based on the roloess; you can input the argument as the bigger the window size, the better fitting the trend line. The last argument is whether to show the plot for the fitted line and offset signal, the 

```matlab
  Winsize = 100;
  Sig1 = Sig1.BgSub(Winsize,'y');
  Back1 = Back1.BgSub(1000,'N');
```
![BgSubCat1](https://user-images.githubusercontent.com/100134089/224961351-3191704a-ea74-41f8-b866-184037d515aa.svg)

## Sampling some peaks based on the height threshold 
This step is based on the blank signal to generate the height threshold and perform sampling follow by this code:

```matlab
Sig1 = Sig1.FlipFindPeak(Back1.Test_signal_offset,'Red','Y',"Offset"); 
```
This method accepts four arguments, the first is the de-trended blank signal, and the next is whether the NIE signal is oxidative or reductive. In 'Red' or 'Oxi', the third argument is to show the plot and show the found peaks on the offset or original trend signal, the output as below:
![FindPeaksCAT](https://user-images.githubusercontent.com/100134089/224968188-8bab673a-283b-4c9e-990f-4a18a26fb868.svg)

## Generating the training sets and Clustering 

Some parameters are extracted to distinguish the types of spikes and noise spikes, performed by the K-means; the method is as 

```matlab
Sig1 = Sig1.GeRawTrainSet('Y','Y');
```
the first argument is based on the elbow method, and the second is by silhouette score, input 'y' or 'n' for the corresponding methods 
#### Elbow method:
![ElbowCAT](https://user-images.githubusercontent.com/100134089/224988058-98917542-1926-47d3-aefd-f1bb2a75ce32.svg)


#### Silhouette Score:
![SC_CAT](https://user-images.githubusercontent.com/100134089/224986902-edcbbe18-0032-42e8-a6cf-07d80d376b41.svg)


## Templates generation
The next step is to generate the representative templates based on the K-means results; the spikes belonging to the same cluster are averaged, obtain the raw templates; the following command makes the method:
```matlab
Sig1 = Sig1.KmeansGeRawSigTem(3,'Y',"Y");
```
it will show the raw templates and bar chart for the sum of different template
![TemRawCAT](https://user-images.githubusercontent.com/100134089/224986136-49ab4ef3-6b7f-46ba-88df-ae843e7fe360.svg)
![BarCAT](https://user-images.githubusercontent.com/100134089/224989387-41b8a0a0-6ba9-401f-b0e7-9ffc79e0acf9.svg)

## Templates Regulation
There are some templates which are templates for the noisy spikes, and some templates need to regulate two sides to avoid involving too much background trend; the first black one is removed by the method output as the last argument works as whether to plot 
```matlab
Sig1 = Sig1.RawtemplatesReguFunc([2,3],'Y');
```
![Regutem](https://user-images.githubusercontent.com/100134089/224990028-70f0515b-8e6c-48ab-9137-b300eacfd2d7.svg)

## Templates Matching 
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
