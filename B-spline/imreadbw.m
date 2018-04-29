function I = imreadbw(file) %对读入的图片如果是RGB的转化为灰度图
% IMREADBW  Reads an image as gray-scale
%   I=IMREADBW(FILE) reads the image FILE and converts the result to a
%   gray scale image (with DOUBLE storage class anr range normalized
%   in [0,1]).

I=im2double(imread(file)) ;
%如果图像img是double型的，d=img;如果图像是logical或single型图像，d=double(img);如果图像是uint8型，d=
%double(img)/255;如果图像是uint16型，d=double(img)/65535;

if(size(I,3) > 1)
  I = rgb2gray( I ) ; %将一副真彩色图像转换成灰度图像
end
%判断图像x是否为RGB图像，如果是，则将图像X转化为灰度图像。
%size(x,3)==3;% 判断是否为RGB图像，是的话转化为灰度图像。
%如果x为RGB图像的话
%i1(:,:,1);% 表示R色
%i1(:,:,2);% 表示G色
%i1(:,:,3);% 表示B色


