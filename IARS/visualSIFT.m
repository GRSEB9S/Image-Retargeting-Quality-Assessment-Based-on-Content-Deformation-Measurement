close all;  
  
% ģ��figure 5  
% im = zeros(401,401,3);  
% im(:,:,:) = 0;   
% im(2:200, 2:200, 2) = 255;   
% im(202:400, 202:400, 2) = 255;   
% im(2:200, 202:400, :) = 255;   
% im(202:400, 2:200, :) = 255;   
% im=im2double(im);  
  
 im = imread('E:\����\datebase\gray\Artroom\ArtRoom.png');  
figure;imshow(im);  
  
cellsize=3;  
gridspacing=1;  
winSize = 3;  
addpath(fullfile(pwd,'MEX','mexDenseSIFT'));  %pwdʶ��ǰ�ļ���fullfile�����ļ���
addpath(fullfile(pwd,'MEX','mexDiscreteFlow'));  
  
%����������(201,201)sift�ֲ�  
sift = mexDenseSIFT(im,cellsize,gridspacing);  
%sift = imnormalize(double(sift));  
x = 1:128;  
y = squeeze(sift(100, 100, :));  
  
figure,plot(x, y);  
xlabel('The dimension of the descriptor in SIFT flow','FontSize',12);  
ylabel('Distribution','FontSize',12);  
set(gca, 'FontSize', 12);  
%����õ�sift�ľ�ֵ�ͷ���  
M1 = mean(y)  
V1 = var(double(y))  
  
% %��������(201,201)sift�ֲ�  
im_noise =imnoise(im, 'gaussian', 0.02);  
sift_noise = mexDenseSIFT(im_noise,cellsize,gridspacing);  
x = 1:128;  
y = squeeze(sift_noise(201, 201, :));  
figure,plot(x, y);  
figure, imshow(im_noise);  
%����õ�sift�ľ�ֵ�ͷ���  
M2 = mean(y)  
V2 = var(double(y))  
  
%����������(100,100)sift�ֲ�  
y1 = squeeze(sift(201, 201, :));  
figure,plot(  x, y1, 'r', x, y, 'b', x, y, '--g', 'LineWidth', 2);  
  
xlabel('Dimension Degree of SIFT Flow Descriptor','FontSize',20);  
ylabel('Distribution','FontSize',20);  
set(gca, 'FontSize', 20);  
axis([0, 130, -0.1, 1.1]);  
  
hleg1 = legend('B', 'A', 'C');  
  
%����õ�sift�ľ�ֵ�ͷ���  
M3 = mean(y)  
V3 = var(double(y))  
  
%%��������(100,100)sift�ֲ�  
y = squeeze(sift_noise(100, 100, :));  
figure,plot(x, y);  
%����õ�sift�ľ�ֵ�ͷ���  
M4 = mean(y)  
V4 = var(double(y))  