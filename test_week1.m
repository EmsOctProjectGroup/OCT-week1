%TODO: load function
%%
figure(1)
colormap gray
image(mscancut)

%%
image(mscancut(:,6150:13750))


%%
arr = zeros(512,6150,66);

bscan_length = 16000;
bscan_start = 6150;
bscan_amount = fix(410001/bscan_length); %410001/bscan_length
bscan_end = (bscan_amount-1)*(bscan_length)
arr = zeros(512,bscan_length,bscan_amount);

index = 1;
for count = bscan_start:bscan_length+bscan_start:bscan_end
    arr(:,:,index) = mscancut(:,count:count+(bscan_length-1));
    index = index+1;
end

bscan0 = mscancut(:,1:6150);
bscan1 = arr(:,:,1);
bscan2 = arr(:,:,2);
bscan3 = arr(:,:,3);
bscan4 = arr(:,:,4);
bscan10 = arr(:,:,25);

figure(2)
colormap gray

subplot(231)
image(bscan0)
title('bscan0')

subplot(232)
image(bscan1)
title('bscan1')

subplot(233)
image(bscan2)
title('bscan2')

subplot(234)
image(bscan3)
title('bscan3')

subplot(235)
image(bscan3)
title('bscan4')

subplot(236)
image(bscan10)
title('bscan10')
%%

bscan1_medfil = medfilt2(bscan1,[2 2])
image(bscan1_medfil)

%%

sob = fspecial('sobel')
bscan1_sob = imfilter(bscan1_medfil,sob)
image(bscan1_sob)
%image(edge(bscan1,'Prewitt'))

bscan1_double = double(bscan1)
scaling = @(x)(x-min(x(:)))/(max(x(:))-min(x(:)))
bscan1_scaleMatrix = scaling(bscan1_double)
bscan1_scaled = (bscan1_scaleMatrix.^10).*200
image(bscan1_scaled)
bscan1_scaled_medfilt = medfilt2(bscan1_scaled, [5 5])
image(bscan1_scaled_medfilt)


%image(edge(bscan1_scaleMatrix,'canny'))

image(bscan2)




