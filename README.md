
# NIE data-driven template matching 

It's a simple toolbox for using the data-driven template matching in nano-impact electrochemistry signal analysis




## Authors

- [@ziwzh166](https://github.com/ziwzh166)
- [@arun-naha](https://github.com/arun-naha)
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
The background suntraction is made by making the smoothdata based on the roloess, you can input the argument as, the bigger the window size the better fitting the trend line, the last argument is wehter or not to show the plot for the fitted line and offset signal

```matlab
  Winsize = 100;
  Sig1 = Sig1.BgSub(Winsize,'y');
  Back1 = Back1.BgSub(1000,'N');
```
![BgSubCat1](https://user-images.githubusercontent.com/100134089/224961351-3191704a-ea74-41f8-b866-184037d515aa.svg)
