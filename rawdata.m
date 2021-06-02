%% Daten Laden
FileName = 'phantom1_2_2raw.dat';
FolderName = 'C:/Users/Bene/Desktop/';
File = fullfile(FolderName, FileName);
fileID = fopen(File)
A = fread(fileID, [1024,410001], 'uint16');

% Matrix A mit 540 multiplizieren
A = A*540;

%% Offset laden
load('C:/Users/Bene/Desktop/offset_chirp.mat', 'Offset');
load('C:/Users/Bene/Desktop/offset_chirp.mat', 'Chirp');

%% Dunkelstrom Kompensieren
B = A-Offset;

%% Plot matrix A und matrix B, um Offset reduziert
subplot(2,2,1)
plot(A(:,1));
subplot(2,2,2)
plot(B(:,1));

%% Entfernung des DC-Terms. Bilde Mittelwert des B-scans und subtrahiere das Ergebnis mit der jeweiligen Spalte
bscan_length=1024
for i = 1:bscan_length
    Bscan_mean(i, 1) = mean(B(i, :));
    B(i, 1) = B(i, 1) - Bscan_mean(i,1);
end

%% Plot ohne DC-Term
plot(B(:,1));

%% Interpolation
interpolation = interp1(1:bscan_length,B(:,1),Chirp);
figure('name','Chirp')
plot(interpolation)

%% Multiplikation mit Fensterfunktion
hw = hann(bscan_length);
multwithhw = interpolation.*hw;
multwithhw(1,1)=0;
figure('name','Hann-Fenster')
plot(multwithhw)

%% Fast Fourier Transformation
fourier = (abs(fft(multwithhw)));
figure('name','A-Scan')
compr = 20*log10(fourier)
plot(compr)
