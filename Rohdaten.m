FileName = 'phantom1_2_2raw.dat';
FolderName = 'C:/Users/haina/Downloads/';
File = fullfile(FolderName, FileName);
fileID = fopen(File);
A=fread(fileID,[1024,410001],'uint16');
%%
A=A*540;
%%
colormap gray;
image(A(:,1000));