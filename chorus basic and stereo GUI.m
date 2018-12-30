classdef chorusBasicStereo_Sonawane_s1889125_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        PlaychoruswaveButton        matlab.ui.control.Button
        ChorusLabel                 matlab.ui.control.Label
        LFO1Label                   matlab.ui.control.Label
        LFO2Label                   matlab.ui.control.Label
        g1StrengthKnobLabel         matlab.ui.control.Label
        g1StrengthKnob              matlab.ui.control.Knob
        D1DelayKnobLabel            matlab.ui.control.Label
        D1DelayKnob                 matlab.ui.control.Knob
        M01DepthKnobLabel           matlab.ui.control.Label
        M01DepthKnob                matlab.ui.control.Knob
        g2StrengthKnobLabel         matlab.ui.control.Label
        g2StrengthKnob              matlab.ui.control.Knob
        D2DelayKnobLabel            matlab.ui.control.Label
        D2DelayKnob                 matlab.ui.control.Knob
        M02DepthKnobLabel           matlab.ui.control.Label
        M02DepthKnob                matlab.ui.control.Knob
        M03DepthKnobLabel           matlab.ui.control.Label
        M03DepthKnob                matlab.ui.control.Knob
        D3DelayKnobLabel            matlab.ui.control.Label
        D3DelayKnob                 matlab.ui.control.Knob
        g3StrengthKnobLabel         matlab.ui.control.Label
        g3StrengthKnob              matlab.ui.control.Knob
        M04DepthKnobLabel           matlab.ui.control.Label
        M04DepthKnob                matlab.ui.control.Knob
        D4DelayKnobLabel            matlab.ui.control.Label
        D4DelayKnob                 matlab.ui.control.Knob
        g4StrengthKnobLabel         matlab.ui.control.Label
        g4StrengthKnob              matlab.ui.control.Knob
        LFO3Label                   matlab.ui.control.Label
        LFO4Label                   matlab.ui.control.Label
        msLabel                     matlab.ui.control.Label
        msLabel_2                   matlab.ui.control.Label
        msLabel_3                   matlab.ui.control.Label
        msLabel_4                   matlab.ui.control.Label
        msLabel_5                   matlab.ui.control.Label
        msLabel_7                   matlab.ui.control.Label
        msLabel_8                   matlab.ui.control.Label
        msLabel_9                   matlab.ui.control.Label
        PlaystereoChoruswaveButton  matlab.ui.control.Button
        LeftChannelMonoLabel        matlab.ui.control.Label
        RightchannelLabel           matlab.ui.control.Label
        LFOFrequency2KnobLabel      matlab.ui.control.Label
        LFOFrequency2Knob           matlab.ui.control.Knob
        LFOFrequency3KnobLabel      matlab.ui.control.Label
        LFOFrequency3Knob           matlab.ui.control.Knob
        LFOFrequency4KnobLabel      matlab.ui.control.Label
        LFOFrequency4Knob           matlab.ui.control.Knob
        LFOFrequency1KnobLabel      matlab.ui.control.Label
        LFOFrequency1Knob           matlab.ui.control.Knob
        HzLabel                     matlab.ui.control.Label
        HzLabel_2                   matlab.ui.control.Label
        HzLabel_3                   matlab.ui.control.Label
        HzLabel_4                   matlab.ui.control.Label
        PlayoriginalwaveButton      matlab.ui.control.Button
    end

    methods (Access = private)

        % Button pushed function: PlaychoruswaveButton
        function PlaychoruswaveButtonPushed(app, event)
            [x,Fs]=audioread('Cath_Very_Short2.wav');
            
            M01value= app.M01DepthKnob.Value;
           M02value=app.M02DepthKnob.Value;
           
           D1value=app.D1DelayKnob.Value;
           D2value=app.D2DelayKnob.Value;
           
           g1value=app.g1StrengthKnob.Value;
           g2value=app.g2StrengthKnob.Value;
           
           f1value=app.LFOFrequency1Knob.Value;
           f2value=app.LFOFrequency2Knob.Value;
           
%Checking for stereo file
if size(x,2)==2                            %detects if stereo or mono
    x = sum(x,2)/size(x,2);                %convert to mono
end
%Length of input signal 
LEN=length(x);           
%--------------------------------------------------------------------------------------%
%LFO1 Parameters
%--------------------------------------------------------------------------------------%
M01=round((M01value/1000)*Fs);    %effect depth for LFO1
D1=round((D1value/1000)*Fs);     %swing range for LFO1
f1=f1value;                   %frequency  for LFO1
g1=g1value;                 %effect strength for LFO1
%--------------------------------------------------------------------------------------%
%LFO2 parameters
%--------------------------------------------------------------------------------------%
M02=round((M02value/1000)*Fs);     %effect depth for LFO1
D2=round((D2value/1000)*Fs);      %swing range for LFO1
f2=f2value;                  %frequency  for LFO1
g2=g2value;                  %effect strength for LFO1
%--------------------------------------------------------------------------------------%
Ts=1/Fs;                %Sample period
%--------------------------------------------------------------------------------------%
%Error Check

if D1>M01 
  error('Depth cant be greater than Delay .Depth =< delay');
  return;
end

if D2>M02
  error('Depth cant be greater than Delay ,Depth =< delay');
  return;
end
%--------------------------------------------------------------------------------------%
%Order of Interpolation
N=2;
%set of values for alpha
Q=6;
%Use of lagrange interpolation coefficient table
coefficient=intertab_sonawane_s1889125(N,Q);
%Initialise range of alpha for the coefficient so as to compare the row
RowOfTable=zeros(1,Q);
 for i=1:Q
   RowOfTable(i)=(-((Q/2) -(i-1))/Q);    %value of alpha between -1/2 to 1/2    
 end
 %Initialise vectors for interpolated values for x(n-M1)
PolyInterpolate1=zeros(length(x),1);
%Initialise vectors for interpolated values for x(n-M2)
PolyInterpolate2=zeros(length(x),1);
%Initialise output signal y
y=zeros(LEN,1);
%INITIALISE VALUES UPTO M01 for output signal
y(1:M01)=x(1:M01);


%--------------------------------------------------------------------------------------%
%Implementation of chorus y(n)= x(n) + g1*x(n-M1(n)) + g2*x(n-M2(n))
%--------------------------------------------------------------------------------------%


for n=M01+1:(LEN)
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
  %Implementation for LFO 1
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%  
    
  %Calculation of M1=D1+M01*sin(2*pi*f1*n*Ts)
   M1(n) =(M01*(sin(2*pi*n*f1*Ts)));
   M1(n)= ((M1(n))+D1);
  %calculation of n-M1
   delta1(n)=n-M1(n);
   
   if delta1(n)-N/2>1
   %Initially value of window for interpolation is consideres as floor of
   %delta
       fdelta1=floor(delta1(n));
         for i=1:N/2
             %Here window means the values to interpolated i.e. if 5.2 is to
             %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
             %Hence window is for the values below floor(delta) and
               %window3 is for the values above floor(delta)
              %Calculation of value for window below floor(delta) 
                window(i)=fdelta1+i;
              %Calculation of value for window above floor(delta)
                window2(i)=fdelta1-i+1;
                window3=flip(window2);
              %Final window to be interpolated for smooth effect
                finalwindow=[window3 window];   %These are polynomial variables 
         end
   %Calculation of alpha
    alpha=delta1(n)-fdelta1-0.5;
    %Initialise index as 0 .This index will be used to compare with the coefficient table
    index1=0;
    difference=abs(RowOfTable-alpha);
    minvalue=min(difference);
    %This is the index of row to be used to get the values of coefficient
    %from the coefficient table of lagrange interpolation
    index1 = find(difference==minvalue);
          for z=1:N
              %Checks if the values are non negative for x(n-M1) and must
              %be within the length of input sigal
              if ((finalwindow>0) & (finalwindow<=LEN))==true
                  %These are the interpolated values for x(n-M1) to be interpolated using the
                  %formula of P(n)=summation of (x(window)*Coefficient)
                 PolyInterpolate1(n)=PolyInterpolate1(n)+(coefficient(index1,z))*x(finalwindow(z));
              end
          end

%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
  %Implementation for LFO 2
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
 %Calculation of M2=D2+M02*sin(2*pi*f2*n*Ts) 
  M2(n) = ( M02*(sin(2*pi*n*f2*Ts)));
  M2(n)= (M2(n))+D2;
  %calculation of n-M1
  delta2(n)=n-M1(n);
   
   if delta2(n)-N/2>0
       %Initially value of window for interpolation is consideres as floor of
       %delta
       fdelta2=floor(delta2(n));
            for i=1:N/2
                %Here window means the values to interpolated i.e. if 5.2 is to
               %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
               %Hence window_2 is for the values below floor(delta) and
               %window3_2 is for the values above floor(delta)
               %Calculation of value for window below floor(delta) 
                     window_2(i)=fdelta2+i;
               %Calculation of value for window above floor(delta)
                     window2_2(i)=fdelta2-i+1;
                     window3_2=flip(window2_2);
                %Final window to be interpolated for smooth effect
                     finalwindow_2=[window3_2 window_2]; %These are polynomial variables 
            end
     %Calculation of alpha
    alpha2=(delta2(n)-fdelta2-0.5);
     %Initialise index as 0 .This index will be used to compare with the coefficient table        
    difference2=abs(RowOfTable-alpha2);
    minvalue2=min(difference2);
    index2 = find(difference2==minvalue2);
            for z=1:N
              %Checks if the values are non negative for x(n-M2)to be interpolated and must
              %be within the length of input sigal
              if ((finalwindow_2>0) & (finalwindow_2<=LEN))==true
                  %These are the interpolated values for x(n-M2) to be interpolated using the
                  %formula of P(n)=summation of (x(window)*Coefficient)
                PolyInterpolate2(n)=PolyInterpolate2(n)+(coefficient(index2,z))*x(finalwindow_2(z));
              end
            end
          
    %Chorus output y=x(n) + x(n-M1)*g1 + x(n-M2)*g2 where n-M1 and n-M2 has the interpolated values. 
       y(n)=x(n)+(PolyInterpolate1(n))*g1+ PolyInterpolate2(n)*g2  ;
        end
    end
end
toc
%Sound of the output signal of chorus
 soundsc(y,Fs)

        end

        % Callback function
        function StereoSwitchValueChanged(app, event)
            value = app.StereoSwitch.Value;
            app.LFO1Label
        end

        % Value changed function: g1StrengthKnob
        function g1StrengthKnobValueChanged(app, event)
            value = app.g1StrengthKnob.Value;
            
        end

        % Value changed function: D1DelayKnob
        function D1DelayKnobValueChanged(app, event)
            value = app.D1DelayKnob.Value;
            
        end

        % Value changed function: M01DepthKnob
        function M01DepthKnobValueChanged(app, event)
            value = app.M01DepthKnob.Value;
            
        end

        % Value changed function: g2StrengthKnob
        function g2StrengthKnobValueChanged(app, event)
            value = app.g2StrengthKnob.Value;
            
        end

        % Value changed function: D2DelayKnob
        function D2DelayKnobValueChanged(app, event)
            value = app.D2DelayKnob.Value;
            
        end

        % Value changed function: M02DepthKnob
        function M02DepthKnobValueChanged(app, event)
            value = app.M02DepthKnob.Value;
            
        end

        % Value changed function: M03DepthKnob
        function M03DepthKnobValueChanged(app, event)
            value = app.M03DepthKnob.Value;
            
        end

        % Value changed function: D3DelayKnob
        function D3DelayKnobValueChanged(app, event)
            value = app.D3DelayKnob.Value;
            
        end

        % Value changed function: g3StrengthKnob
        function g3StrengthKnobValueChanged(app, event)
            value = app.g3StrengthKnob.Value;
            
        end

        % Value changed function: M04DepthKnob
        function M04DepthKnobValueChanged(app, event)
            value = app.M04DepthKnob.Value;
            
        end

        % Value changed function: D4DelayKnob
        function D4DelayKnobValueChanged(app, event)
            value = app.D4DelayKnob.Value;
            
        end

        % Value changed function: g4StrengthKnob
        function g4StrengthKnobValueChanged(app, event)
            value = app.g4StrengthKnob.Value;
            
        end

        % Button pushed function: PlaystereoChoruswaveButton
        function PlaystereoChoruswaveButtonPushed(app, event)
             [x,Fs]=audioread('Cath_Very_Short2.wav');
            
           M01value= app.M01DepthKnob.Value;
           M02value=app.M02DepthKnob.Value;
           M03value=app.M03DepthKnob.Value;
           M04value=app.M04DepthKnob.Value;
       
           
           D1value=app.D1DelayKnob.Value;
           D2value=app.D2DelayKnob.Value;
           D3value=app.D3DelayKnob.Value;
           D4value=app.D4DelayKnob.Value;
           
           g1value=app.g1StrengthKnob.Value;
           g2value=app.g2StrengthKnob.Value;
           g3value=app.g3StrengthKnob.Value;
           g4value=app.g3StrengthKnob.Value;
           
           f1value=app.LFOFrequency1Knob.Value;
           f2value=app.LFOFrequency2Knob.Value;
           f3value=app.LFOFrequency3Knob.Value;
           f4value=app.LFOFrequency4Knob.Value;
           
           %Checking for stereo file
if size(x,2)==2                            %detects if stereo or mono
    x = sum(x,2)/size(x,2);                %convert to mono
end
%Length of input signal 
LEN=length(x); 

%Declare LFO Parameters
%--------------------------------------------------------------------------------------%
%LFO1 Parameters
%--------------------------------------------------------------------------------------%
M01=round((M01value/1000)*Fs);    %effect depth for LFO1
D1=round((D1value/1000)*Fs);     %swing range for LFO1
f1=f1value;                   %frequency  for LFO1
g1=g1value;                 %effect strength for LFO1
%--------------------------------------------------------------------------------------%
%LFO2 Parameters
%--------------------------------------------------------------------------------------%
M02=round((M02value/1000)*Fs);     %effect depth for LFO2
D2=round((D2value/1000)*Fs);      %swing range for LFO2
f2=f2value;                    %frequency  for LFO2
g2=g2value;                  %effect strength for LFO2
%--------------------------------------------------------------------------------------%
%LFO1 Parameters
%--------------------------------------------------------------------------------------%
M03=round((M03value/1000)*Fs);    %effect depth for LFO1
D3=round((D3value/1000)*Fs);     %swing range for LFO1
f3=f3value;                   %frequency  for LFO1
g3=g3value;                 %effect strength for LFO2
%--------------------------------------------------------------------------------------%
%LFO2 Parameters
%--------------------------------------------------------------------------------------%
M04=round((M04value/1000)*Fs);     %effect depth for LFO1
D4=round((D4value/1000)*Fs);      %swing range for LFO1
f4=f4value;                    %frequency  for LFO1
g4=g4value;                 %effect strength for LFO2
%--------------------------------------------------------------------------------------%
Ts=1/Fs;                  %Sample period
%--------------------------------------------------------------------------------------%

%Error Check
%--------------------------------------------------------------------------------------%
%Error check for LFO1 Parameters
%--------------------------------------------------------------------------------------%
if D1>M01 
  error('Depth cant be greater than Delay.  Depth =< Delay');
  return;
end
if (g1<= -1 | g1>=1)==true
    error('Effect strenght cannot be greater than 1 and less than -1');
    return;
end
%--------------------------------------------------------------------------------------%
%Error check for LFO2 Parameters
%--------------------------------------------------------------------------------------%

if D2>M02
  error('Depth cant be greater than Delay .Depth =< delay');
  return;
end
if (g2<= -1 | g2>=1)==true
    error('Effect strenght cannot be greater than 1 and less than -1');
    return;
end

%--------------------------------------------------------------------------------------%
%Error check for LFO3 Parameters
%--------------------------------------------------------------------------------------%

if D3>M03 
  error('Depth cant be greater than Delay.  Depth =< Delay');
  return;
end
if (g3<= -1 | g3>=1)==true
    error('Effect strenght cannot be greater than 1 and less than -1');
    return;
end
%--------------------------------------------------------------------------------------%
%Error check for LFO4 Parameters
%--------------------------------------------------------------------------------------%
if D4>M04
  error('Depth cant be greater than Delay . Depth =< delay');
  return;
end
if (g4<= -1 | g4>=1)==true
    error('Effect strenght cannot be greater than 1 and less than -1');
    return;
end

%--------------------------------------------------------------------------------------%
%Order of Interpolation
N=2;
%set of values for alpha
Q=4;
%Use of lagrange interpolation coefficient table
RowOfTable=zeros(1,Q);
%Use of lagrange interpolation coefficient table
coefficient=intertab_sonawane_s1889125(N,Q);
%Initialise range of alpha for the coefficient so as to compare the row
 for i=1:Q
   RowOfTable(i)=(-((Q/2) -(i-1))/Q);    %value of alpha between -1/2 to 1/2    
 end
%Initialise vectors for interpolated values for x(n-M1)
PolyInterpolate1=zeros(length(x),1);
%Initialise vectors for interpolated values for x(n-M2)
PolyInterpolate2=zeros(length(x),1);
%Initialise vectors for interpolated values for x(n-M3)
PolyInterpolate3=zeros(length(x),1);
%Initialise vectors for interpolated values for x(n-M4)
PolyInterpolate4=zeros(length(x),1);
%Initialise output signal for left channel 
yL=zeros(LEN,1);
%INITIALISE VALUES UPTO M01 for output signal
yL(1:M01)=x(1:M01);
%Initialise output signal for right channel 
yR=zeros(LEN,1);
%INITIALISE VALUES UPTO M02 for output signal
yR(1:M02)=x(1:M02);


%--------------------------------------------------------------------------------------%
%*Implementation of chorus for left channel yL(n)= x(n) + g1*x(n-M1(n)) + g2*x(n-M2(n))
%--------------------------------------------------------------------------------------%

for n=M01+1:(LEN)
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
  %Implementation for LFO 1
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
  
%Calculation of M1=D1+M01*sin(2*pi*f1*n*Ts)
M1(n) =(M01*(sin(2*pi*n*f1*Ts)));
M1(n)= ((M1(n))+D1);
%calculation of n-M1
delta1(n)=n-M1(n);
   
   if delta1(n)-N/2>1
   %Initially value of window for interpolation is consideres as floor of
   %delta
       fdelta1=floor(delta1(n));
           for i=1:N/2
               %Here window means the values to interpolated i.e. if 5.2 is to
             %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
             %Hence window is for the values below floor(delta) and
               %window3 is for the values above floor(delta)
              %Calculation of value for window below floor(delta)
                window(i)=fdelta1+i;
                 %Calculation of value for window above floor(delta)
                window2(i)=fdelta1-i+1;
                window3=flip(window2);
                 %Final window to be interpolated for smooth effect
                finalwindow=[window3 window];   % %These are polynomial variables 
           end
     %Calculation of alpha
      alpha=delta1(n)-fdelta1-0.5;
      %Initialise index as 0 .This index will be used to compare with the coefficient table
      index1=0;
      difference=abs(RowOfTable-alpha);
      minvalue=min(difference);
      index1 = find(difference==minvalue);
       %This is the index of row to be used to get the values of coefficient
    %from the coefficient table of lagrange interpolation
          for z=1:N
              %Checks if the values are non negative for x(n-M1) and must
              %be within the length of input sigal
              if ((finalwindow>0) & (finalwindow<=LEN))==true
                  %These are the interpolated values for x(n-M1) to be interpolated using the
                  %formula of P(n)=summation of (x(window)*Coefficient)
             PolyInterpolate1(n)=PolyInterpolate1(n)+(coefficient(index1,z))*x(finalwindow(z));
              end
          end

%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
     %Implementation for LFO 2
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
%Calculation of M2=D2+M02*sin(2*pi*f2*n*Ts)  
M2(n) = ( M02*(sin(2*pi*n*f2*Ts)));
M2(n)= (M2(n))+D2;
%calculation of n-M1
delta2(n)=n-M1(n);
   
   if delta2(n)-N/2>0
       %Initially value of window for interpolation is consideres as floor of
       %delta
       fdelta2=floor(delta2(n));
            for i=1:N/2
                %Here window means the values to interpolated i.e. if 5.2 is to
               %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
               %Hence window_2 is for the values below floor(delta) and
               %window3_2 is for the values above floor(delta)
               %Calculation of value for window below floor(delta) 
                     window_2(i)=fdelta2+i;
                     %Calculation of value for window above floor(delta)
                     window2_2(i)=fdelta2-i+1;
                     window3_2=flip(window2_2);
                     %Final window to be interpolated for smooth effect
                     finalwindow_2=[window3_2 window_2]; %These are polynomial variables
            end
    %Calculation of alpha
    alpha2=(delta2(n)-fdelta2-0.5);
    %Initialise index as 0 .This index will be used to compare with the coefficient table         
    difference2=abs(RowOfTable-alpha2);
    minvalue2=min(difference2);
    index2 = find(difference2==minvalue2);
             for z=1:N
               %Checks if the values are non negative for x(n-M2)to be interpolated and must
              %be within the length of input sigal
                 if ((finalwindow_2>0) & (finalwindow_2<=LEN))==true
                     %These are the interpolated values for x(n-M2) to be interpolated using the
                     %formula of P(n)=summation of (x(window)*Coefficient)
                    PolyInterpolate2(n)=PolyInterpolate2(n)+(coefficient(index2,z))*x(finalwindow_2(z));
                 end
             end
            %Chorus output for left channel yL=x(n) + x(n-M1)*g1 + x(n-M2)*g2 where n-M1 and n-M2 has the interpolated values.
       yL(n)=x(n)+(PolyInterpolate1(n))*g1+ PolyInterpolate2(n)*g2  ;
        end
    end
  end
%    end
%end
yL=transpose(yL);

%--------------------------------------------------------------------------------------%
%Implementation of chorus for Right channel yR(n)= x(n) + g3*x(n-M3(n)) + g4*x(n-M4(n))
%--------------------------------------------------------------------------------------%


for n=M03+1:(LEN)
    
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
  %Implementation for LFO 3
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%  
%Calculation of M1=D1+M01*sin(2*pi*f1*n*Ts)   
M3(n) =(M03*(cos(2*pi*n*f1*Ts)));
M3(n)= ((M3(n))+D3);
%calculation of n-M3
delta3(n)=n-M3(n);
   
   if delta3(n)-N/2>1
   %Initially value of window for interpolation is consideres as floor of
   %delta3
       fdelta3=floor(delta3(n));
         for i=1:N/2
             %Here window means the values to interpolated i.e. if 5.2 is to
             %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
             %Hence window is for the values below floor(delta) and
               %window3 is for the values above floor(delta)
              %Calculation of value for window below floor(delta) 
                 window_3(i)=fdelta3+i;
                 %Calculation of value for window above floor(delta)
                 window2_3(i)=fdelta3-i+1;
                 window3_3=flip(window2_3);
                 %Final window to be interpolated for smooth effect
                 finalwindow_3=[window3_3 window_3];   %These are polynomial variables
         end
    %Calculation of alpha
      alpha3=delta3(n)-fdelta3-0.5;
    %Initialise index as 0 .This index will be used to compare with the coefficient table
      index3=0;
      difference3=abs(RowOfTable-alpha3);
      minvalue3=min(difference3);
      %This is the index of row to be used to get the values of coefficient
      %from the coefficient table of lagrange interpolation
      index3 = find(difference3==minvalue3);
            for z=1:N
                 %Checks if the values are non negative for x(n-M3) and must
                 %be within the length of input sigal
                    if ((finalwindow_3>0) & (finalwindow_3<=LEN))==true
                        %These are the interpolated values for x(n-M1) to be
                        %interpolated using the  formula of P(n)=summation of (x(window)*Coefficient)
                         PolyInterpolate3(n)=PolyInterpolate3(n)+(coefficient(index1,z))*x(finalwindow_3(z));
                    end
             end

%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
%Implementation for LFO 4
%--------------------------------------------------------------------------------------%
%--------------------------------------------------------------------------------------%
%Calculation of M4=D4+M04*sin(2*pi*f4*n*Ts) 
M4(n) = ( M04*(cos(2*pi*n*f2*Ts)));
M4(n)= (M4(n))+D4;
%calculation of n-M1
delta4(n)=n-M4(n);
   
   if delta4(n)-N/2>0
       %Initially value of window for interpolation is consideres as floor of
       %delta
       fdelta4=floor(delta4(n));
            for i=1:N/2
                %Here window means the values to interpolated i.e. if 5.2 is to
               %be interpolated about N=4 then finalwindow will be  (4,5,6,7)
               %Hence window_2 is for the values below floor(delta) and
               %window3_2 is for the values above floor(delta)
               %Calculation of value for window below floor(delta) 
                     window_4(i)=fdelta4+i;
                     %Calculation of value for window above floor(delta)
                     window2_4(i)=fdelta4-i+1;
                     window3_4=flip(window2_4);
                     %Final window to be interpolated for smooth effect
                     finalwindow_4=[window3_4 window_4];
            end
   %Calculation of alpha
    alpha4=(delta4(n)-fdelta4-0.5);
    %Initialise index as 0 .This index will be used to compare with the coefficient table 
    difference4=abs(RowOfTable-alpha4);
    minvalue4=min(difference4);
    index4 = find(difference4==minvalue4);
            for z=1:N
               %Checks if the values are non negative for x(n-M2)to be interpolated and must
              %be within the length of input sigal
                if ((finalwindow_4>0) & (finalwindow_4<=LEN))==true
                  %These are the interpolated values for x(n-M2) to be interpolated using the
                  %formula of P(n)=summation of (x(window)*Coefficient)
                     PolyInterpolate4(n)=PolyInterpolate4(n)+(coefficient(index2,z))*x(finalwindow_4(z));
                end
            end
 %Chorus output for right channel yR=x(n) + x(n-M3)*g3 + x(n-M4)*g4 where n-M3 and n-M4 has the interpolated values.
       yR(n)=x(n)+(PolyInterpolate3(n))*g3+ PolyInterpolate4(n)*g4  ;
              end
   end
        end
yR=transpose(yR);
%Stereo chorus output signal y
y=[yL;yR];
toc
soundsc(y,Fs)
           
           
           
           
           
           
        end

        % Value changed function: LFOFrequency2Knob
        function LFOFrequency2KnobValueChanged(app, event)
            value = app.LFOFrequency2Knob.Value;
            
        end

        % Value changed function: LFOFrequency3Knob
        function LFOFrequency3KnobValueChanged(app, event)
            value = app.LFOFrequency3Knob.Value;
            
        end

        % Value changed function: LFOFrequency4Knob
        function LFOFrequency4KnobValueChanged(app, event)
            value = app.LFOFrequency4Knob.Value;
            
        end

        % Value changed function: LFOFrequency1Knob
        function LFOFrequency1KnobValueChanged(app, event)
            value = app.LFOFrequency1Knob.Value;
            
        end

        % Button pushed function: PlayoriginalwaveButton
        function PlayoriginalwaveButtonPushed(app, event)
            [x,Fs]=audioread('Cath_Very_Short2.wav');
            soundsc(x,Fs);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 970 808];
            app.UIFigure.Name = 'UI Figure';

            % Create PlaychoruswaveButton
            app.PlaychoruswaveButton = uibutton(app.UIFigure, 'push');
            app.PlaychoruswaveButton.ButtonPushedFcn = createCallbackFcn(app, @PlaychoruswaveButtonPushed, true);
            app.PlaychoruswaveButton.FontName = 'Arial Black';
            app.PlaychoruswaveButton.Position = [788 371 131 46];
            app.PlaychoruswaveButton.Text = 'Play chorus wave';

            % Create ChorusLabel
            app.ChorusLabel = uilabel(app.UIFigure);
            app.ChorusLabel.FontName = 'Brush Script MT';
            app.ChorusLabel.FontSize = 66;
            app.ChorusLabel.FontWeight = 'bold';
            app.ChorusLabel.FontAngle = 'italic';
            app.ChorusLabel.Position = [334 736 155 87];
            app.ChorusLabel.Text = 'Chorus';

            % Create LFO1Label
            app.LFO1Label = uilabel(app.UIFigure);
            app.LFO1Label.FontName = 'Damascus';
            app.LFO1Label.FontAngle = 'italic';
            app.LFO1Label.Position = [109 692 38 22];
            app.LFO1Label.Text = 'LFO 1';

            % Create LFO2Label
            app.LFO2Label = uilabel(app.UIFigure);
            app.LFO2Label.FontAngle = 'italic';
            app.LFO2Label.Position = [306 692 38 22];
            app.LFO2Label.Text = 'LFO 2';

            % Create g1StrengthKnobLabel
            app.g1StrengthKnobLabel = uilabel(app.UIFigure);
            app.g1StrengthKnobLabel.HorizontalAlignment = 'center';
            app.g1StrengthKnobLabel.Position = [95 205 68 22];
            app.g1StrengthKnobLabel.Text = 'g1 Strength';

            % Create g1StrengthKnob
            app.g1StrengthKnob = uiknob(app.UIFigure, 'continuous');
            app.g1StrengthKnob.Limits = [0.1 0.9];
            app.g1StrengthKnob.MajorTicks = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            app.g1StrengthKnob.ValueChangedFcn = createCallbackFcn(app, @g1StrengthKnobValueChanged, true);
            app.g1StrengthKnob.MinorTicks = [];
            app.g1StrengthKnob.Position = [96 261 65 65];
            app.g1StrengthKnob.Value = 0.2;

            % Create D1DelayKnobLabel
            app.D1DelayKnobLabel = uilabel(app.UIFigure);
            app.D1DelayKnobLabel.HorizontalAlignment = 'center';
            app.D1DelayKnobLabel.Position = [97 371 54 22];
            app.D1DelayKnobLabel.Text = 'D1 Delay';

            % Create D1DelayKnob
            app.D1DelayKnob = uiknob(app.UIFigure, 'continuous');
            app.D1DelayKnob.Limits = [1 9];
            app.D1DelayKnob.ValueChangedFcn = createCallbackFcn(app, @D1DelayKnobValueChanged, true);
            app.D1DelayKnob.MinorTicks = [];
            app.D1DelayKnob.Position = [93 427 60 60];
            app.D1DelayKnob.Value = 9;

            % Create M01DepthKnobLabel
            app.M01DepthKnobLabel = uilabel(app.UIFigure);
            app.M01DepthKnobLabel.HorizontalAlignment = 'center';
            app.M01DepthKnobLabel.Position = [94 535 65 22];
            app.M01DepthKnobLabel.Text = 'M01 Depth';

            % Create M01DepthKnob
            app.M01DepthKnob = uiknob(app.UIFigure, 'continuous');
            app.M01DepthKnob.Limits = [11 39];
            app.M01DepthKnob.MajorTicks = [11 16 21 26 31 36 39];
            app.M01DepthKnob.ValueChangedFcn = createCallbackFcn(app, @M01DepthKnobValueChanged, true);
            app.M01DepthKnob.MinorTicks = [7 12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 38];
            app.M01DepthKnob.Position = [96 591 60 60];
            app.M01DepthKnob.Value = 28;

            % Create g2StrengthKnobLabel
            app.g2StrengthKnobLabel = uilabel(app.UIFigure);
            app.g2StrengthKnobLabel.HorizontalAlignment = 'center';
            app.g2StrengthKnobLabel.Position = [283 207 68 22];
            app.g2StrengthKnobLabel.Text = 'g2 Strength';

            % Create g2StrengthKnob
            app.g2StrengthKnob = uiknob(app.UIFigure, 'continuous');
            app.g2StrengthKnob.Limits = [0.1 0.9];
            app.g2StrengthKnob.MajorTicks = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            app.g2StrengthKnob.ValueChangedFcn = createCallbackFcn(app, @g2StrengthKnobValueChanged, true);
            app.g2StrengthKnob.MinorTicks = [];
            app.g2StrengthKnob.Position = [286 263 60 60];
            app.g2StrengthKnob.Value = 0.1;

            % Create D2DelayKnobLabel
            app.D2DelayKnobLabel = uilabel(app.UIFigure);
            app.D2DelayKnobLabel.HorizontalAlignment = 'center';
            app.D2DelayKnobLabel.Position = [288 371 54 22];
            app.D2DelayKnobLabel.Text = 'D2 Delay';

            % Create D2DelayKnob
            app.D2DelayKnob = uiknob(app.UIFigure, 'continuous');
            app.D2DelayKnob.Limits = [1 9];
            app.D2DelayKnob.MajorTicks = [1 2 3 4 5 6 7 8 9];
            app.D2DelayKnob.ValueChangedFcn = createCallbackFcn(app, @D2DelayKnobValueChanged, true);
            app.D2DelayKnob.Position = [284 427 60 60];
            app.D2DelayKnob.Value = 4;

            % Create M02DepthKnobLabel
            app.M02DepthKnobLabel = uilabel(app.UIFigure);
            app.M02DepthKnobLabel.HorizontalAlignment = 'center';
            app.M02DepthKnobLabel.Position = [282 535 65 22];
            app.M02DepthKnobLabel.Text = 'M02 Depth';

            % Create M02DepthKnob
            app.M02DepthKnob = uiknob(app.UIFigure, 'continuous');
            app.M02DepthKnob.Limits = [11 39];
            app.M02DepthKnob.MajorTicks = [11 16 21 26 31 36 39];
            app.M02DepthKnob.ValueChangedFcn = createCallbackFcn(app, @M02DepthKnobValueChanged, true);
            app.M02DepthKnob.MinorTicks = [7 12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 38];
            app.M02DepthKnob.Position = [284 591 60 60];
            app.M02DepthKnob.Value = 17;

            % Create M03DepthKnobLabel
            app.M03DepthKnobLabel = uilabel(app.UIFigure);
            app.M03DepthKnobLabel.HorizontalAlignment = 'center';
            app.M03DepthKnobLabel.Position = [469 535 65 22];
            app.M03DepthKnobLabel.Text = 'M03 Depth';

            % Create M03DepthKnob
            app.M03DepthKnob = uiknob(app.UIFigure, 'continuous');
            app.M03DepthKnob.Limits = [11 39];
            app.M03DepthKnob.MajorTicks = [11 16 21 26 31 34 39];
            app.M03DepthKnob.ValueChangedFcn = createCallbackFcn(app, @M03DepthKnobValueChanged, true);
            app.M03DepthKnob.MinorTicks = [7 12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 38];
            app.M03DepthKnob.Position = [471 591 60 60];
            app.M03DepthKnob.Value = 11;

            % Create D3DelayKnobLabel
            app.D3DelayKnobLabel = uilabel(app.UIFigure);
            app.D3DelayKnobLabel.HorizontalAlignment = 'center';
            app.D3DelayKnobLabel.Position = [484 371 54 22];
            app.D3DelayKnobLabel.Text = 'D3 Delay';

            % Create D3DelayKnob
            app.D3DelayKnob = uiknob(app.UIFigure, 'continuous');
            app.D3DelayKnob.Limits = [1 9];
            app.D3DelayKnob.MajorTicks = [1 2 3 4 5 6 7 8 9];
            app.D3DelayKnob.ValueChangedFcn = createCallbackFcn(app, @D3DelayKnobValueChanged, true);
            app.D3DelayKnob.Position = [480 427 60 60];
            app.D3DelayKnob.Value = 1;

            % Create g3StrengthKnobLabel
            app.g3StrengthKnobLabel = uilabel(app.UIFigure);
            app.g3StrengthKnobLabel.HorizontalAlignment = 'center';
            app.g3StrengthKnobLabel.Position = [485 210 68 22];
            app.g3StrengthKnobLabel.Text = 'g3 Strength';

            % Create g3StrengthKnob
            app.g3StrengthKnob = uiknob(app.UIFigure, 'continuous');
            app.g3StrengthKnob.Limits = [0.1 0.9];
            app.g3StrengthKnob.MajorTicks = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            app.g3StrengthKnob.ValueChangedFcn = createCallbackFcn(app, @g3StrengthKnobValueChanged, true);
            app.g3StrengthKnob.MinorTicks = [];
            app.g3StrengthKnob.Position = [488 266 60 60];
            app.g3StrengthKnob.Value = 0.4;

            % Create M04DepthKnobLabel
            app.M04DepthKnobLabel = uilabel(app.UIFigure);
            app.M04DepthKnobLabel.HorizontalAlignment = 'center';
            app.M04DepthKnobLabel.Position = [660 535 65 22];
            app.M04DepthKnobLabel.Text = 'M04 Depth';

            % Create M04DepthKnob
            app.M04DepthKnob = uiknob(app.UIFigure, 'continuous');
            app.M04DepthKnob.Limits = [11 39];
            app.M04DepthKnob.MajorTicks = [11 16 21 26 31 34 39];
            app.M04DepthKnob.ValueChangedFcn = createCallbackFcn(app, @M04DepthKnobValueChanged, true);
            app.M04DepthKnob.MinorTicks = [7 12 13 14 15 17 18 19 20 22 23 24 25 27 28 29 30 32 33 34 35 38];
            app.M04DepthKnob.Position = [662 591 60 60];
            app.M04DepthKnob.Value = 11;

            % Create D4DelayKnobLabel
            app.D4DelayKnobLabel = uilabel(app.UIFigure);
            app.D4DelayKnobLabel.HorizontalAlignment = 'center';
            app.D4DelayKnobLabel.Position = [665 371 54 22];
            app.D4DelayKnobLabel.Text = 'D4 Delay';

            % Create D4DelayKnob
            app.D4DelayKnob = uiknob(app.UIFigure, 'continuous');
            app.D4DelayKnob.Limits = [1 9];
            app.D4DelayKnob.MajorTicks = [1 2 3 4 5 6 7 8 9];
            app.D4DelayKnob.ValueChangedFcn = createCallbackFcn(app, @D4DelayKnobValueChanged, true);
            app.D4DelayKnob.Position = [662 427 60 60];
            app.D4DelayKnob.Value = 1;

            % Create g4StrengthKnobLabel
            app.g4StrengthKnobLabel = uilabel(app.UIFigure);
            app.g4StrengthKnobLabel.HorizontalAlignment = 'center';
            app.g4StrengthKnobLabel.FontAngle = 'italic';
            app.g4StrengthKnobLabel.Position = [668 210 68 22];
            app.g4StrengthKnobLabel.Text = 'g4 Strength';

            % Create g4StrengthKnob
            app.g4StrengthKnob = uiknob(app.UIFigure, 'continuous');
            app.g4StrengthKnob.Limits = [0.1 0.9];
            app.g4StrengthKnob.MajorTicks = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
            app.g4StrengthKnob.ValueChangedFcn = createCallbackFcn(app, @g4StrengthKnobValueChanged, true);
            app.g4StrengthKnob.MinorTicks = [];
            app.g4StrengthKnob.FontAngle = 'italic';
            app.g4StrengthKnob.Position = [672 266 60 60];
            app.g4StrengthKnob.Value = 0.1;

            % Create LFO3Label
            app.LFO3Label = uilabel(app.UIFigure);
            app.LFO3Label.Position = [482 692 38 22];
            app.LFO3Label.Text = 'LFO 3';

            % Create LFO4Label
            app.LFO4Label = uilabel(app.UIFigure);
            app.LFO4Label.Position = [688 692 38 22];
            app.LFO4Label.Text = 'LFO 4';

            % Create msLabel
            app.msLabel = uilabel(app.UIFigure);
            app.msLabel.Position = [115 556 25 22];
            app.msLabel.Text = 'ms';

            % Create msLabel_2
            app.msLabel_2 = uilabel(app.UIFigure);
            app.msLabel_2.Position = [112 392 25 22];
            app.msLabel_2.Text = 'ms';

            % Create msLabel_3
            app.msLabel_3 = uilabel(app.UIFigure);
            app.msLabel_3.Position = [679 392 25 22];
            app.msLabel_3.Text = 'ms';

            % Create msLabel_4
            app.msLabel_4 = uilabel(app.UIFigure);
            app.msLabel_4.Position = [495 556 25 22];
            app.msLabel_4.Text = 'ms';

            % Create msLabel_5
            app.msLabel_5 = uilabel(app.UIFigure);
            app.msLabel_5.Position = [503 392 25 22];
            app.msLabel_5.Text = 'ms';

            % Create msLabel_7
            app.msLabel_7 = uilabel(app.UIFigure);
            app.msLabel_7.Position = [303 556 25 22];
            app.msLabel_7.Text = 'ms';

            % Create msLabel_8
            app.msLabel_8 = uilabel(app.UIFigure);
            app.msLabel_8.Position = [299 392 25 22];
            app.msLabel_8.Text = 'ms';

            % Create msLabel_9
            app.msLabel_9 = uilabel(app.UIFigure);
            app.msLabel_9.Position = [682 556 25 22];
            app.msLabel_9.Text = 'ms';

            % Create PlaystereoChoruswaveButton
            app.PlaystereoChoruswaveButton = uibutton(app.UIFigure, 'push');
            app.PlaystereoChoruswaveButton.ButtonPushedFcn = createCallbackFcn(app, @PlaystereoChoruswaveButtonPushed, true);
            app.PlaystereoChoruswaveButton.FontName = 'Arial Black';
            app.PlaystereoChoruswaveButton.Position = [784 141 178 35];
            app.PlaystereoChoruswaveButton.Text = 'Play stereo Chorus wave';

            % Create LeftChannelMonoLabel
            app.LeftChannelMonoLabel = uilabel(app.UIFigure);
            app.LeftChannelMonoLabel.FontName = 'Arial Rounded MT Bold';
            app.LeftChannelMonoLabel.FontWeight = 'bold';
            app.LeftChannelMonoLabel.Position = [168 702 138 61];
            app.LeftChannelMonoLabel.Text = 'Left Channel / Mono';

            % Create RightchannelLabel
            app.RightchannelLabel = uilabel(app.UIFigure);
            app.RightchannelLabel.FontName = 'Arial Rounded MT Bold';
            app.RightchannelLabel.FontWeight = 'bold';
            app.RightchannelLabel.Position = [560 721 86 22];
            app.RightchannelLabel.Text = 'Right channel';

            % Create LFOFrequency2KnobLabel
            app.LFOFrequency2KnobLabel = uilabel(app.UIFigure);
            app.LFOFrequency2KnobLabel.HorizontalAlignment = 'center';
            app.LFOFrequency2KnobLabel.Position = [269 32 98 22];
            app.LFOFrequency2KnobLabel.Text = 'LFO Frequency 2';

            % Create LFOFrequency2Knob
            app.LFOFrequency2Knob = uiknob(app.UIFigure, 'continuous');
            app.LFOFrequency2Knob.Limits = [0.1 4.9];
            app.LFOFrequency2Knob.ValueChangedFcn = createCallbackFcn(app, @LFOFrequency2KnobValueChanged, true);
            app.LFOFrequency2Knob.Position = [287 88 60 60];
            app.LFOFrequency2Knob.Value = 0.8;

            % Create LFOFrequency3KnobLabel
            app.LFOFrequency3KnobLabel = uilabel(app.UIFigure);
            app.LFOFrequency3KnobLabel.HorizontalAlignment = 'center';
            app.LFOFrequency3KnobLabel.Position = [462 32 98 22];
            app.LFOFrequency3KnobLabel.Text = 'LFO Frequency 3';

            % Create LFOFrequency3Knob
            app.LFOFrequency3Knob = uiknob(app.UIFigure, 'continuous');
            app.LFOFrequency3Knob.Limits = [0.1 4.9];
            app.LFOFrequency3Knob.ValueChangedFcn = createCallbackFcn(app, @LFOFrequency3KnobValueChanged, true);
            app.LFOFrequency3Knob.Position = [480 88 60 60];
            app.LFOFrequency3Knob.Value = 0.4;

            % Create LFOFrequency4KnobLabel
            app.LFOFrequency4KnobLabel = uilabel(app.UIFigure);
            app.LFOFrequency4KnobLabel.HorizontalAlignment = 'center';
            app.LFOFrequency4KnobLabel.Position = [668 32 98 22];
            app.LFOFrequency4KnobLabel.Text = 'LFO Frequency 4';

            % Create LFOFrequency4Knob
            app.LFOFrequency4Knob = uiknob(app.UIFigure, 'continuous');
            app.LFOFrequency4Knob.Limits = [0.1 4.9];
            app.LFOFrequency4Knob.ValueChangedFcn = createCallbackFcn(app, @LFOFrequency4KnobValueChanged, true);
            app.LFOFrequency4Knob.Position = [686 88 60 60];
            app.LFOFrequency4Knob.Value = 0.8;

            % Create LFOFrequency1KnobLabel
            app.LFOFrequency1KnobLabel = uilabel(app.UIFigure);
            app.LFOFrequency1KnobLabel.HorizontalAlignment = 'center';
            app.LFOFrequency1KnobLabel.Position = [79 32 98 22];
            app.LFOFrequency1KnobLabel.Text = 'LFO Frequency 1';

            % Create LFOFrequency1Knob
            app.LFOFrequency1Knob = uiknob(app.UIFigure, 'continuous');
            app.LFOFrequency1Knob.Limits = [0.1 4.9];
            app.LFOFrequency1Knob.ValueChangedFcn = createCallbackFcn(app, @LFOFrequency1KnobValueChanged, true);
            app.LFOFrequency1Knob.Position = [97 88 60 60];
            app.LFOFrequency1Knob.Value = 0.4;

            % Create HzLabel
            app.HzLabel = uilabel(app.UIFigure);
            app.HzLabel.Position = [115 53 25 22];
            app.HzLabel.Text = 'Hz';

            % Create HzLabel_2
            app.HzLabel_2 = uilabel(app.UIFigure);
            app.HzLabel_2.Position = [305 53 25 22];
            app.HzLabel_2.Text = 'Hz';

            % Create HzLabel_3
            app.HzLabel_3 = uilabel(app.UIFigure);
            app.HzLabel_3.Position = [498 53 25 22];
            app.HzLabel_3.Text = 'Hz';

            % Create HzLabel_4
            app.HzLabel_4 = uilabel(app.UIFigure);
            app.HzLabel_4.Position = [704 53 25 22];
            app.HzLabel_4.Text = 'Hz';

            % Create PlayoriginalwaveButton
            app.PlayoriginalwaveButton = uibutton(app.UIFigure, 'push');
            app.PlayoriginalwaveButton.ButtonPushedFcn = createCallbackFcn(app, @PlayoriginalwaveButtonPushed, true);
            app.PlayoriginalwaveButton.FontName = 'Arial Black';
            app.PlayoriginalwaveButton.Position = [784 591 135 41];
            app.PlayoriginalwaveButton.Text = 'Play original wave';
        end
    end

    methods (Access = public)

        % Construct app
        function app = chorusBasicStereo_Sonawane_s1889125_GUI

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
