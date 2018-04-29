% Example of nonrigid registration using steepest gradient optimizerʹ������ݶ��Ż����ķǸ�����׼��ʾ��
% Example is written by D.Kroon University of Twente (October 2008)
% clean
clear all; close all; 

% Read two greyscale images of Lena
%  I1=im2double(imread('lenag1.png'));
%  I2=im2double(imread('lenag2.png'));
I1=im2double(imread('E:\HYT\NRID\ABdata_N\ours_3\ours_3_0.75_multiop.jpg'));
I2=im2double(imread('E:\HYT\NRID\ABdata_N\ours_3\ours_3.jpg'));

% Type of registration error used see registration_error.m
options.type='sd';  %FIR����Ӧ�˲���
options.centralgrad=false;   % Use fast forward instead of central error gradient
Spacing=[8 8];  % b-spline grid spacing in x and y direction�����С��32*32
[O_trans]=make_init_grid(Spacing,size(I1));  %�� % Make the Initial b-spline registration grid������ʼ��׼����,����Ǿ�������
I1=double(I1); I2=double(I2); O_trans=double(O_trans);   % Convert all values to type double

% Smooth both images for faster registration
I1s=imfilter(I1,fspecial('gaussian'));
I2s=imfilter(I2,fspecial('gaussian'));

% Optimizer parameters 
%��Display','iter',��ʾÿһ���ĵ�����'GradObj','on'��ʾ�ݶ�ֵ��'MaxIter'����������������'DiffMinChange'
% DiffMaxChange �C �������޲���ݶȵ����仯��DiffMinChange - �������޲���ݶȵ���С�仯��
optim=struct('Display','iter','GradObj','on','MaxIter',20,'DiffMinChange',0.1,'DiffMaxChange',1);

% Reshape O_trans from a matrix to a vector.�������ת��Ϊ����
sizes=size(O_trans); O_trans=O_trans(:);

% Start the b-spline nonrigid registration optimizer���ǵ�������
% fminsd�ҵ����������ĺ����ľֲ���Сֵ��ʹ����Ľ����½��Ż���
% X = fminsd��FUN��X0����X0����ʼ�������Բ��Ҿֲ���Сֵ
%  X�ĺ���FUN�� FUN��������X������һ������
%  ��X�����Ƶĺ���ֵF.X0�����Ǳ�������������󡣴��Ż������ٶȿ���ͨ���ṩ�������X�����ݶȡ�
%  fminsd ��ȡ���淽��ʽ�ľֲ���Сֵ ����bspline_registration_gradient   ��Դ˾��Ż� help fminsd��ϸ
O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim); %��

% Reshape O_trans from a vector to a matrix
O_trans=reshape(O_trans,sizes);

% Transform the input image with the found optimal grid.����bspline_transform
Icor=bspline_transform(O_trans,I1,Spacing); %��

% Make a (transformed) grid image
Igrid=make_grid_image(Spacing,size(I1));  %��
% figure,
% imshow(Igrid);hold on;
[Igrid,Tx,Ty,Tmap]=bspline_transform(O_trans,Igrid,Spacing);  %��jw[I,Tx,Ty,Tmap,Tz

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
