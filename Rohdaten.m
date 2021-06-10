FileName = 'phantom1_2_2raw.dat';
FolderName = 'C:/Users/haina/Downloads/';
File = fullfile(FolderName, FileName);
fileID = fopen(File);
A=fread(fileID,[1024,inf],'uint16');
A=A*540;

%% Offset
B_off=A-Offset;

%% DC-Term entfernen
DC =  B_off - mean(B_off,2);

%% Interpolation
interpolation=interp1(DC,Chirp);

%%
subplot(2,1,1);
plot(DC(:,1));
subplot(2,1,2);
plot(interpolation(:,1));

%% Hann
hw = hann(1024);
multwithhw = interpolation.*hw;
multwithhw(1,:)=0;
%%
plot(multwithhw(:,187));
%% Fourier
fourier=abs(fft(multwithhw));
compr=20*log10(fourier);

%%
compr_mult=fix(compr*2.46);
%%
figure(1);
colormap gray(350);
image(compr_mult(1:512,:));
% figure(2);
% colormap gray(350);
% image(mscancut);
