%%
% loading mscancut
FileName = 'phantom1_2_2.mat';
FolderName = '~/Desktop/Fr_Gr_5';
File = fullfile(FolderName, FileName);

load(File);

%%
% print loaded data as picture - grayscale
figure(1)
colormap gray
image(mscancut)

%%
image(mscancut(:,6150:13750))

%%
arr = zeros(512,6150,66);
bscan_start = 6150;
%%

bscan_length = 7800;
bscan_start = 6180;
bscan_amount = fix((410001-bscan_start)/bscan_length); %410001/bscan_length
bscan_end = (bscan_amount-1)*(bscan_length)

arr = zeros(512,bscan_length,bscan_amount);

index = 1;
for count = bscan_start:bscan_length:bscan_start+(bscan_amount-1)*bscan_length
    arr(:,:,index) = mscancut(:,count:count+(bscan_length-1));
    index = index+1;
end
%%
bscan_amount

%%
bscan1 = arr(:,:,1);

image(bscan1)


%%
for x = 1:bscan_amount
    subplot(8,7,x);
    image(arr(:,:,x));
end

%%

matmean = double(mean(bscan1,'all'))

%%
colormap gray
for i = 0:200
    for j = 1:512:512*7650
        index = i+j;
        bscan1(index) = matmean;
        bscan1(221+j) = matmean;
        bscan1(222+j) = matmean;
    end
end
image(bscan1)

%%

figure(2)
colormap gray
image(bscan1)

%%
figure(3)
colormap gray
bscan_gauss = conv2(bscan1,fspecial('gaussian', [5 5], 1));
image(bscan_gauss)

%%
figure(4)
colormap gray
bscan_medfilt = medfilt2(bscan_gauss,[5 5]);
image(bscan_medfilt)

%%
figure(5)
colormap gray
bscan_medfilt_imguided = imguidedfilter(bscan_medfilt);
image(bscan_medfilt_imguided)

%%

matmean = double(mean(bscan_medfilt_imguided,'all'))

%%
figure(6)
colormap gray
bscan_medfilt_imguided_contrast = zeros(512,7650);
for x = 1:512
    for y = 1:7650
        if bscan_medfilt_imguided(x,y) > matmean+14
            bscan_medfilt_imguided_contrast(x,y) = bscan_medfilt_imguided(x,y)+150;
        else
            bscan_medfilt_imguided_contrast(x,y) = bscan_medfilt_imguided(x,y)-100;
        end
    end
end
image(bscan_medfilt_imguided_contrast)

%%

colormap gray

for x = 1:512*7650-1
    bscan_medfilt_imguided_contrast(x) = bscan_medfilt_imguided_contrast(x+1)-bscan_medfilt_imguided_contrast(x);
end

image(bscan_medfilt_imguided_contrast)


%%
colormap gray
bscan_wiener = wiener2(bscan_medfilt_imguided_contrast,[5 5]);
image(bscan_wiener)

%%
colormap gray
sob = fspecial('sobel')
bscan1_sob = imfilter(bscan_medfilt_imguided_contrast,sob);
image(bscan1_sob)

%%
colormap(gray(2))
bscan_edge = edge(bscan1_sob,'Prewitt');
image(bscan_edge)
%%
figure(7)
colormap gray
bscan_medfilt_contrast_imguided_stdfilt = stdfilt(bscan1_sob);
image(bscan_medfilt_contrast_imguided_stdfilt)

%%
im2 = bscan_medfilt_imguided_contrast;

[nrows, ncols] = size(im2);

increment = 2*pi/7650;
startAngle = 0;

rho = repmat([1:nrows]',1,ncols);

theta = repmat([startAngle:increment:startAngle + increment*(ncols-1)],nrows,1);

[x,y] = pol2cart(theta, rho);

[zz, xx, yy] = ffgrid(x,y,im2,1,1);

Z = griddata(x,y,im2,xx,yy');
%%
colormap gray
image(Z)

%%
colormap gray

nans = isnan(Z);

nans(1)

for x = 1:1024*1025
    if nans(x)
        Z(x) = matmean;
    end
end

image(Z)

    

%%
function pcimg=imgpolarcoord(img,radius,angle)
% IMGPOLARCOORD converts a given image from cartesian coordinates to polar
% coordinates.
%
% Input:
%        img  : bidimensional image.
%      radius : radius length (# of pixels to be considered).
%      angle  : # of angles to be considered for decomposition.
%
% Output:
%       pcimg : polar coordinate image.
%
% Usage Example:
%  im=double(imread('cameraman.tif'));
%  fim=fft2(im);
%  pcimg=iapolarcoord(im);
%  fpcimg=iapolarcoord(fim);
%  figure; subplot(2,2,1); imagesc(im); colormap gray; axis image;
%  title('Input image');  subplot(2,2,2);
%  imagesc(log(abs(fftshift(fim)+1)));  colormap gray; axis image;
%  title('FFT');subplot(2,2,3); imagesc(pcimg); colormap gray; axis image;
%  title('Polar Input image');  subplot(2,2,4);
%  imagesc(log(abs(fpcimg)+1));  colormap gray; axis image;
%  title('Polar FFT');
%
% Notes:
%  The software is provided "as is", without warranty of any kind.
%  Javier Montoya would like to thank prof. Juan Carlos Gutierrez for his
%  support and suggestions, while studying polar-coordinates.
%  Authors: Juan Carlos Gutierrez & Javier Montoya.
   if nargin < 1
      error('Please specify an image!');
   end
   
   img         = double(img);
   [rows,cols] = size(img);
   cy          = round(rows/2);
   cx          = round(cols/2);
   
   if exist('radius','var') == 0
      radius = min(round(rows/2),round(cols/2))-1;
   end
   
   if exist('angle','var') == 0
      angle = 360;
   end
  
   pcimg = [];
   i     = 1;
   
   for r=0:radius
      j = 1;
      for a=0:2*pi/angle:2*pi-2*pi/angle
         pcimg(i,j) = img(cy+round(r*sin(a)),cx+round(r*cos(a)));
         j = j + 1;
      end
      i = i + 1;
   end
end

function [zzgrid, xxvec, yyvec] = ffgrid(x, y, z, dx, dy, xyz0);
% [zgrid, xvec, yvec] = ffgrid(x, y, z, dx, dy, [x0 x1 y0 y1 z0 z1 n0]);
%
%  Fast 'n' Furious data gridding.
%
%  Input:  Unevenly spaced data in vectors x, y, and z of equal length.
%          If dx and dy are omitted, the default is 75 bins in each direction.
%
%          If dx or dy is negative, then the variable is taken as the number of
%          bins rather than a grid resolution. If dx is complex, the imaginary
%          part is taken as the value with which empty grid points are padded.
%          The default padding is 10% below minimum value in grid.
%
%          The vector containing the limits can be padded with NaNs if only
%          certain limits are desired, e g if x1 and z0 are wanted:
%
%            ffgrid(x, y, z, [.5 nan nan nan 45])
%
%          The last parameter, n0, removes outliers from the data set by 
%          ignoring grid points with n0 or less observations. When n0 is 
%          negative it is treated as the percentage of the total number of 
%          data points.
%
%  Output: The gridded matrix ZGRID together with vectors XVEC and YVEC. 
%          If no output arguments are given, FFGRID will plot the gridded
%          function with the prescribed axes using PCOLOR.
%
%  Requires bin.m. Tested under MatLab 4.2, 5.0, and 5.1.
%
%  See also bin.m for further details about dx, x0, etc, and density.m.
%

% 28.7.97, Oyvind.Breivik@gfi.uib.no.
%
% Oyvind Breivik
% Department of Geophysics
% University of Bergen
% NORWAY

DX = -75; % default value

x = x(:);

y = y(:);

z = z(:);

xyz = NaN*ones(1,7);

if (nargin < 6)
 xyz0 = min(x);
end

if (nargin == 4 & length(dx) > 1)
 xyz0 = dx;
 dx = DX;
end
if (nargin < 4)
 dx = DX;
end
if (real(dx) == 0)
 dx = DX + dx;
end
pad = imag(dx);
dx = real(dx);

if (nargin == 5 & length(dy) > 1)
  xyz0 = dy;
  dy = dx;
end
if (nargin < 5)
  dy = dx;
end

nxyz = length(xyz0);
xyz(1:nxyz) = xyz0;

if (isnan(xyz(7)))
  xyz(7) = 0;
end
if (isnan(xyz(6)))
  xyz(6) = max(z);
end
if (isnan(xyz(5)))
  xyz(5) = min(z);
end
if (isnan(xyz(4)))
  xyz(4) = max(y);
end
if (isnan(xyz(3)))
  xyz(3) = min(y);
end
if (isnan(xyz(2)))
  xyz(2) = max(x);
end
if (isnan(xyz(1)))
  xyz(1) = min(x);
end
x0 = xyz(1); x1 = xyz(2); y0 = xyz(3); y1 = xyz(4); z0 = xyz(5); z1 = xyz(6);
n0 = xyz(7);
  
if (dx < 0)
 dx = (x1 - x0)/abs(dx);
end
if (dy < 0)
 dy = (y1 - y0)/abs(dy);
end

ix = bin(x, dx, x0, x1);
iy = bin(y, dy, y0, y1); % bin data in (x,y)-space

xvec = x0:dx:x1;
yvec = y0:dy:y1;

nx = length(xvec);
ny = length(yvec);

inx = (ix >= 1) & (ix <= nx);
iny = (iy >= 1) & (iy <= ny);
inz = (z >= z0) & (z <= z1);
in = (inx & iny & inz);
ix = ix(in); iy = iy(in); z = z(in);
N = length(ix); % how many datapoints are left now?

ngrid = zeros(nx, ny); % no of obs per grid cell
zgrid = ngrid; % z-coordinate

for i = 1:N
 zgrid(ix(i), iy(i)) = zgrid(ix(i), iy(i)) + z(i);
 ngrid(ix(i), iy(i)) = ngrid(ix(i), iy(i)) + 1;
end

% Remove outliers
if (n0 >= 0)
  zgrid(ngrid <= n0) = 0;
  ngrid(ngrid <= n0) = 0;
else
  n0 = -n0;
  zgrid(ngrid <= n0*N) = 0;
  ngrid(ngrid <= n0*N) = 0;
end

N = sum(sum(ngrid)); % how many datapoints are left now?

Nil = (ngrid == 0);
ngrid(Nil) = 1; % now we don't divide by zero

zgrid = zgrid./ngrid; % Make average of values for each point

if (~pad)
  zmax = max(max(zgrid(~Nil)));
  zmin = min(min(zgrid(~Nil)));
  pad = zmin - (zmax - zmin)/10; % Adjust padding value to values in grid
end

zgrid(Nil) = pad; % Empty grid points are set to default value
rho = nnz(~Nil)/length(Nil(:));

zgrid = zgrid'; % Get in shape

if (nargout == 0) % no output, then plot
 pcolor(xvec, yvec, zgrid)
 colorbar
 colormap jet
 xlabel(inputname(1))
 ylabel(inputname(2))
 zstr=inputname(3);
 dum = size(zgrid');
 if (~isempty(zstr)) % all this vital information ...
  str = sprintf('Color scale: %s, %d data points, grid: %dx%d, density: %4.2f', ...
  inputname(3), N, dum(1)-1, dum(2)-1, rho);
  title(str);
 else
   str = sprintf('%d data points, grid: %dx%d, density: %4.2f', ...
         N, dum(1)-1, dum(2)-1, rho);
   title(str);
 end
end

if (nargout > 0)
 zzgrid = zgrid;
end

if (nargout > 1)
 xxvec = xvec;
 yyvec = yvec;
end
end
function [i, nnbins] = bin(x, dx, x0, x1);
%
% [i, nbins] = bin(x, dx, x0, x1);
%
% Returns the vector of indices, starting from 1, 
% corresponding to the chosen bin size, dx, 
% start x0 and end x1. If x1 is omitted, x1 = max(x) - dx/2. 
% If x0 is omitted, x0 = min(x) + dx/2. If dx is omitted, the data 
% are divided into 10 classes. Note that outliers are not removed.
%
% Tested under MatLab 4.2, 5.0, and 5.1.
%

% 17.1.97, Oyvind.Breivik@gfi.uib.no.
%
% Oyvind Breivik
% Department of Geophysics
% University of Bergen
% NORWAY

N = 10; % Default is 10 classes

if nargin < 2
 dx = (max(x) - min(x))/N;
end
if nargin < 3
 x0 = min(x) + dx/2;
end
if nargin < 4
 x1 = max(x) - dx/2;
end
nbins = round((x1 - x0)/dx) + 1;
i = round((x - x0)/dx) + 1;
%in = (i >= 1) & (i <= nbins); % Indices are within range [1, nbins]. 
%i = i(in);

if nargout > 1
 nnbins = nbins;
end
end

function J=regiongrowing(I,x,y,reg_maxdist)
% This function performs "region growing" in an image from a specified
% seedpoint (x,y)
%
% J = regiongrowing(I,x,y,t) 
% 
% I : input image 
% J : logical output image of region
% x,y : the position of the seedpoint (if not given uses function getpts)
% t : maximum intensity distance (defaults to 0.2)
%
% The region is iteratively grown by comparing all unallocated neighbouring pixels to the region. 
% The difference between a pixel's intensity value and the region's mean, 
% is used as a measure of similarity. The pixel with the smallest difference 
% measured this way is allocated to the respective region. 
% This process stops when the intensity difference between region mean and
% new pixel become larger than a certain treshold (t)
%
% Example:
%
% I = im2double(imread('medtest.png'));
% x=198; y=359;
% J = regiongrowing(I,x,y,0.2); 
% figure, imshow(I+J);
%
% Author: D. Kroon, University of Twente

if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end

J = zeros(size(I)); % Output 
Isizes = size(I); % Dimensions of input image

reg_mean = I(x,y); % The mean of the segmented region
reg_size = 1; % Number of pixels in region

% Free memory to store neighbours of the (segmented) region
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 

pixdist=0; % Distance of the region newest pixel to the regio mean

% Neighbor locations (footprint)
neigb=[-1 0; 1 0; 0 -1;0 1];

% Start regiogrowing until distance between regio and posible new pixels become
% higher than a certain treshold
while(pixdist<reg_maxdist&&reg_size<numel(I))

    % Add new neighbors pixels
    for j=1:4,
        % Calculate the neighbour coordinate
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        
        % Check if neighbour is inside or outside the image
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        
        % Add neighbor if inside and not already part of the segmented area
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

    % Add a new block of free memory
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
    % Add pixel with intensity nearest to the mean of the region, to the region
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
    % Calculate the new mean of the region
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = neg_list(index,1); y = neg_list(index,2);
    
    % Remove the pixel from the neighbour (check) list
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end

% Return the segmented area as logical matrix
J=J>1;
end
