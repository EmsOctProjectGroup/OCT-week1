% FileName = 'phantom1_2_2raw.dat';
% FolderName = 'C:/Users/haina/Downloads/';
% File = fullfile(FolderName, FileName);

fileID = fopen(File);
A=fread(fileID,[1024,inf],'uint16');
A=A*540;

%% Offset
A=A-Offset;

%% DC-Term entfernen
A =  A - mean(A,2);

%% Interpolation
A=interp1(Chirp,A,0:1023);

%% Hann
hw = hann(1024);
A = A.*hw;
%% Fourier
A=abs(fft(A));
A=20*log(A);

%%
A=fix(A*1.069);

function processed = processraw(file,length)
    fileID = fopen(file);
    A=fread(fileID,[length,inf],'uint16');
    A=A*540;

%% Offset
	A=A-Offset;

%% DC-Term entfernen
    A =  A - mean(A,2);

%% Interpolation
    A=interp1(Chirp,A,0:length-1);

%% Hann
    hw = hann(length);
    A = A.*hw;
%% Fourier
    A=abs(fft(A));
    A=20*log(A);

%%
    A=fix(A*1.069);
    processed=A;
end
