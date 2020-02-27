%reads data from the file named filename, and returns sampled data, x, and a sample rate for that data, Fs.
[x, Fs]=audioread('x1.wav');    

%Enter a positive number for up shifting       (for example pitchnote=1  ,increase 1 semitone)
%Enter a negative number for down shifting     (for example pitchnote=-1 ,decrease 1 semitone)
pitchnote=-24;

pitch=pitchnote*(Fs/2);
if pitch<=0;   
    ZerosNumber=-pitch;
    DeleteNumber=0;
else
    DeleteNumber=pitch;
    ZerosNumber=0;
end
     
%Time calculation
N=length(x); 
t=(0:1/Fs:(N-1)/Fs);

%Graphic illustration of the original signal
subplot(211) 
plot(t,x),title('Time Domain Speach Signal'),xlabel('Time (t) '),ylabel('Amplitude')

%Fourier Transform of the signal
SignalFft=fft(x,N); 

%Down Shifting
AddZeros = [zeros(ZerosNumber,1);SignalFft];         %it adds zeros to the signal
SignalFft=AddZeros;

%Up Shifting
DeleteRow=SignalFft;
DeleteRow(1:DeleteNumber,:)=[];                      %it deletes to the signal's rows (1:number of rows that it delete,:)=[];
SignalFft=DeleteRow;


%Frequency calculation of the FFT signal
Np=length(SignalFft);
f2=Fs/Np.*(0:Np-1); 

%Graphic illustration of the FFT signal
subplot(212) 
plot(f2,SignalFft),title('Frequency Spectrum'),xlabel('Frequency (Hz)'),ylabel('Ampltide');

%Inverse Fourier Transform of the signal
SignalIfft=ifft(SignalFft);

%Output of the pitch shifted signal
Outputx=real(SignalIfft);

%sound(Outputx,Fs);                                   %listen to the signal
filename='x2.wav';
audiowrite(filename,Outputx,Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
factor=length(x)/(length(x)-pitch);
%factor=2;
[signal,Fs]=audioread('x2.wav');


winsize=1024;
fftsize=winsize;
window=hann(winsize);
shop=winsize/4;
ahop=floor(shop/factor);

SignalLen=length(signal);
num_win = floor((SignalLen-winsize)/ahop);

OutLen=(num_win-1)*shop+winsize;
Out=zeros(OutLen,1);

TWOPI=2*pi;
CENTERFREQ = [0:fftsize-1]*TWOPI*ahop/fftsize; %Piece of equation showed in q
first=1;
PosIn=1;
PosOut=1;

for win_count=1:num_win
   if first==1
      framed=signal(1:winsize).*window;
      X=fft(fftshift(framed),fftsize);
      Mag=abs(X);
      Pha=angle(X);
      PhaSy=Pha;            
      first=0;
   else
      framed = signal(PosIn:PosIn+winsize-1) .* window; %framed and windowed / current position analysis
      X=fft(fftshift(framed),fftsize);  %Apply FFT whith circular shift
      Mag=abs(X); % Get the Magnitude
      Pha=angle(X); %Get the Phase
      phaseDiff = Pha - old_pha; %Difference between the current and previous phase 
      phaseDiff = phaseDiff - CENTERFREQ'; %Expected phase (unwrapped phase)
      dphi = phaseDiff - TWOPI * round(phaseDiff /TWOPI); %principal argument, MAP phase to +/- pi
      freq = (CENTERFREQ + dphi') /ahop; %true frequency
      PhaSy = old_PhaSy + shop*freq'; %Phase synthesis

   end

   Y = Mag .* ( cos(PhaSy) + sqrt(-1) *(sin(PhaSy)) ); %Resynthesis
   y_out = fftshift(real(ifft(Y,fftsize))).*window; %back to time domain 

   Out(PosOut:PosOut+winsize-1)= Out(PosOut:PosOut+winsize-1) + y_out; %Overlap and add

   old_pha=Pha;
   old_PhaSy=PhaSy;
   PosIn = PosIn + ahop;
   PosOut=PosOut+shop;
end

%Play
sound(Out,Fs);