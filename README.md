



# NIE data-driven template matching 

It's a simple toolbox for using the data-driven template matching in nano-impact electrochemistry signal analysis




## Authors

- [@ziwzh166](https://github.com/ziwzh166)
- [@arun-naha](https://github.com/arun-naha)



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

In the folder /test, you can find two .m files, SignalAnalysis is the toolbox for performing the analysis, and the Display.m is the live editor file which you can see the results step by step. and the other four .txt files are the NIE signals and blank signals we mentioned in the paper.

In the fisrt section are created to imprt the data and creaste new classes
which is given by, the first name string in the folder is imported into the class.

```matlab
    Sig1 = SignalAnalysis(names_Str(1));
```

## Resampling 

As we state in the paper, for the uneven sampled signal, it's integrated in the preprocess method

```matlab
  Sig1 = Sig1.Preprocess("PA","N","Y","Y");
```

the first argument is made for change the current unit, this simple toolbox will convert all scale of current into pA, the second arugment is wether or not to shift the time from 0, can be input as 'yes' or 'no'. then the third arugment is the resample is wether performed or not, the last one is if to show the plot with origianl and resampled signal. 
It will be the plot as below: 
![Resampling](https://user-images.githubusercontent.com/100134089/224944121-91084fef-a1f0-4e92-b4d2-900f48043e30.svg)

## Denoise 
As we state in the paper, we denoise the signal by the 
