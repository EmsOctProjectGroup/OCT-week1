%%
colorbar
colormap gray(255)
bscan1 = mscancut(:,1:6150);

bscan1_float = double(bscan1);

%[a,b] nach [c, d]
% U = (U -a) / (b-a)(d-c) +c
scaling = @(x) (x - min(x(:))) / (max(x(:)) - min(x(:)));
scaled_img = scaling(bscan1_float).*255;

filtered_img = medfilt2(scaled_img, [5 5]);
%filtered_img = medfilt2(filtered_img, [5 5]);
subplot 141
imagesc(filtered_img)

subplot 142
imagesc(bscan1)

size(bscan1)

artifact = zeros(512, 6150);

for i = 220:225
    artifact(i,:) = 120;
end

sob = fspecial('sobel');
one_sob = imfilter(filtered_img, sob);
two_sob = imfilter(one_sob, sob');

imgHrVer = max(one_sob, two_sob);

subplot 143
imagesc(one_sob)

subplot 144
imagesc(imgHrVer)
%%

subplot 143
image(artifact)

subplot 144
image(filtered_img - artifact)

%%
