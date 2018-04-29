function I = imreadbw(file) %�Զ����ͼƬ�����RGB��ת��Ϊ�Ҷ�ͼ
% IMREADBW  Reads an image as gray-scale
%   I=IMREADBW(FILE) reads the image FILE and converts the result to a
%   gray scale image (with DOUBLE storage class anr range normalized
%   in [0,1]).

I=im2double(imread(file)) ;
%���ͼ��img��double�͵ģ�d=img;���ͼ����logical��single��ͼ��d=double(img);���ͼ����uint8�ͣ�d=
%double(img)/255;���ͼ����uint16�ͣ�d=double(img)/65535;

if(size(I,3) > 1)
  I = rgb2gray( I ) ; %��һ�����ɫͼ��ת���ɻҶ�ͼ��
end
%�ж�ͼ��x�Ƿ�ΪRGBͼ������ǣ���ͼ��Xת��Ϊ�Ҷ�ͼ��
%size(x,3)==3;% �ж��Ƿ�ΪRGBͼ���ǵĻ�ת��Ϊ�Ҷ�ͼ��
%���xΪRGBͼ��Ļ�
%i1(:,:,1);% ��ʾRɫ
%i1(:,:,2);% ��ʾGɫ
%i1(:,:,3);% ��ʾBɫ


