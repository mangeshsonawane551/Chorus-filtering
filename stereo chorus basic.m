%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Program Details: An implementation of  Stereo chorus application of the varying  
%delay lengths, M1[n] and M2[n], as well as M3[n] and M4[n] that uses Lagrange interpola-
%tion to ensure a smooth effect. 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

clear all;
clc;
tic
%Reads input audio file and its sample rate
[x,Fs]=audioread('Cath_Very_Short2.wav');
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
M01=round(0.030*Fs);    %effect depth for LFO1
D1=round(0.009*Fs);     %swing range for LFO1
f1=0.2;                   %frequency  for LFO1
g1=0.4;                 %effect strength for LFO1
%--------------------------------------------------------------------------------------%
%LFO2 Parameters
%--------------------------------------------------------------------------------------%
M02=round(0.017*Fs);     %effect depth for LFO2
D2=round(0.004*Fs);      %swing range for LFO2
f2=1;                    %frequency  for LFO2
g2=0.1;                  %effect strength for LFO2
%--------------------------------------------------------------------------------------%
%LFO1 Parameters
%--------------------------------------------------------------------------------------%
M03=round(0.030*Fs);    %effect depth for LFO1
D3=round(0.009*Fs);     %swing range for LFO1
f3=0.2;                   %frequency  for LFO1
g3=0.4;                 %effect strength for LFO2
%--------------------------------------------------------------------------------------%
%LFO2 Parameters
%--------------------------------------------------------------------------------------%
M04=round(0.017*Fs);     %effect depth for LFO1
D4=round(0.004*Fs);      %swing range for LFO1
f4=1;                    %frequency  for LFO1
g4=0.1;                 %effect strength for LFO2
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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
