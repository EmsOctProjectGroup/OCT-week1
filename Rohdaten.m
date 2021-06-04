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
for count = 1:410001
    interpolation(:,count)=interp1(1:1024,DC(:,count),Chirp);
end
%%
plot(interpolation(:,1));
%% Hann
hw = hann(1024);
multwithhw = interpolation.*hw;
multwithhw(1,:)=0;
figure('name','Hann-Fenster')
%%
plot(multwithhw(:,1256));
%% Fourier
for count = 1:410001
    fourier(:,count) = (abs(fft(multwithhw(:,count))));
    compr(:,count) = 20*log10(fourier(:,count));
end
%%
compr_mult=fix(compr*2.46);
%%
figure(1);
colormap gray(350);
image(compr_mult(1:512,:));
figure(2);
colormap gray(350);
image(mscancut);
