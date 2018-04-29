% Example of nonrigid registration using steepest gradient optimizer使用最陡的梯度优化器的非刚性配准的示例
% Example is written by D.Kroon University of Twente (October 2008)
% clean
clear all; close all; 

% Read two greyscale images of Lena
%  I1=im2double(imread('lenag1.png'));
%  I2=im2double(imread('lenag2.png'));
I1=im2double(imread('E:\HYT\NRID\ABdata_N\ours_3\ours_3_0.75_multiop.jpg'));
I2=im2double(imread('E:\HYT\NRID\ABdata_N\ours_3\ours_3.jpg'));

% Type of registration error used see registration_error.m
options.type='sd';  %FIR自适应滤波器
options.centralgrad=false;   % Use fast forward instead of central error gradient
Spacing=[8 8];  % b-spline grid spacing in x and y direction网格大小是32*32
[O_trans]=make_init_grid(Spacing,size(I1));  %调 % Make the Initial b-spline registration grid创建初始配准网格,输出是均匀网格
I1=double(I1); I2=double(I2); O_trans=double(O_trans);   % Convert all values to type double

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial('gaussian'));
I2s=imfilter(I2,fspecial('gaussian'));

% Optimizer parameters 
%‘Display','iter',显示每一步的迭代；'GradObj','on'显示梯度值；'MaxIter'最大允许迭代次数；'DiffMinChange'
% DiffMaxChange – 变量有限差分梯度的最大变化。DiffMinChange - 变量有限差分梯度的最小变化。
optim=struct('Display','iter','GradObj','on','MaxIter',20,'DiffMinChange',0.1,'DiffMaxChange',1);

% Reshape O_trans from a matrix to a vector.网格矩阵转换为向量
sizes=size(O_trans); O_trans=O_trans(:);

% Start the b-spline nonrigid registration optimizer正是迭代过程
% fminsd找到几个变量的函数的局部最小值，使用最陡的渐变下降优化。
% X = fminsd（FUN，X0）在X0处开始，并尝试查找局部最小值
%  X的函数FUN。 FUN接受输入X并返回一个标量
%  在X处估计的函数值F.X0可以是标量，向量或矩阵。此优化器的速度可以通过提供来提高在X处的梯度。
%  fminsd 求取后面方程式的局部最小值 调用bspline_registration_gradient   针对此句优化 help fminsd详细
O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim); %调

% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);

% Transform the input image with the found optimal grid.调用bspline_transform
Icor=bspline_transform(O_trans,I1,Spacing); %调

% Make a (transformed) grid image
Igrid=make_grid_image(Spacing,size(I1));  %调
% figure,
% imshow(Igrid);hold on;
[Igrid,Tx,Ty,Tmap]=bspline_transform(O_trans,Igrid,Spacing);  %调jw[I,Tx,Ty,Tmap,Tz

% Show the registration results
figure,
subplot(2,2,1), imshow(I1); title('input image 1');
subplot(2,2,2), imshow(I2); title('input image 2');
subplot(2,2,3), imshow(Icor); title('transformed image 1');
subplot(2,2,4), imshow(Igrid); title('grid');

figure,
% imshow(I1);hold on;
imshow(Igrid); title('grid');
figure,imshow(Tmap);
