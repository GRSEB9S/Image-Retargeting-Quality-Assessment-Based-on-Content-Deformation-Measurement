close all;  
  
% 模拟figure 5  
% im = zeros(401,401,3);  
% im(:,:,:) = 0;   
% im(2:200, 2:200, 2) = 255;   
% im(202:400, 202:400, 2) = 255;   
% im(2:200, 202:400, :) = 255;   
% im(202:400, 2:200, :) = 255;   
% im=im2double(im);  
  
 im = imread('E:\代码\datebase\gray\Artroom\ArtRoom.png');  
figure;imshow(im);  
  
cellsize=3;  
gridspacing=1;  
winSize = 3;  
addpath(fullfile(pwd,'MEX','mexDenseSIFT'));  %pwd识别当前文件，fullfile创建文件名
addpath(fullfile(pwd,'MEX','mexDiscreteFlow'));  
  
%不加噪声，(201,201)sift分布  
sift = mexDenseSIFT(im,cellsize,gridspacing);  
%sift = imnormalize(double(sift));  
x = 1:128;  
y = squeeze(sift(100, 100, :));  
  
figure,plot(x, y);  
xlabel('The dimension of the descriptor in SIFT flow','FontSize',12);  
ylabel('Distribution','FontSize',12);  
set(gca, 'FontSize', 12);  
%计算该点sift的均值和方差  
M1 = mean(y)  
V1 = var(double(y))  
  
% %加噪声，(201,201)sift分布  
im_noise =imnoise(im, 'gaussian', 0.02);  
sift_noise = mexDenseSIFT(im_noise,cellsize,gridspacing);  
x = 1:128;  
y = squeeze(sift_noise(201, 201, :));  
figure,plot(x, y);  
figure, imshow(im_noise);  
%计算该点sift的均值和方差  
M2 = mean(y)  
V2 = var(double(y))  
  
%不加噪声，(100,100)sift分布  
y1 = squeeze(sift(201, 201, :));  
figure,plot(  x, y1, 'r', x, y, 'b', x, y, '--g', 'LineWidth', 2);  
  
xlabel('Dimension Degree of SIFT Flow Descriptor','FontSize',20);  
ylabel('Distribution','FontSize',20);  
set(gca, 'FontSize', 20);  
axis([0, 130, -0.1, 1.1]);  
  
hleg1 = legend('B', 'A', 'C');  
  
%计算该点sift的均值和方差  
M3 = mean(y)  
V3 = var(double(y))  
  
%%加噪声，(100,100)sift分布  
y = squeeze(sift_noise(100, 100, :));  
figure,plot(x, y);  
%计算该点sift的均值和方差  
M4 = mean(y)  
V4 = var(double(y))  