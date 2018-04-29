% Example of nonrigid registration using steepest gradient optimizerʹ������ݶ��Ż����ķǸ�����׼��ʾ��
% Example is written by D.Kroon University of Twente (October 2008)
clear all; close all; 

%% path and other initial information
PATH_ROOT = 'D:\HYT\gray\'; % the path direct to the MIT dataset
load([PATH_ROOT 'subjData-ref_37.mat'])
subj_data = subjData.data;
SET_NUM = 37; % 37 set images
C1 = 1e-6;   
PATH_NAME = cell(SET_NUM,1);        %37��һ�еĿվ���
All_ratio = zeros(SET_NUM,1);       %37��һ�е������
for set_num = 1:SET_NUM            %ѭ��37��    
    foo = subjData.datasetNames{set_num}; % ��{}����ͼ�������
    foo_loc = strfind(foo,'_0.');       %Ѱ��foo�д��С����е�λ��
    PATH_NAME{set_num} =  foo(1:foo_loc-1);   %·��������ʾ��'_0'ǰһλ
    if(strfind(foo,'_0.75'))           %��ʾ���ű�
        All_ratio(set_num) = 75;     
    elseif(strfind(foo,'_0.50'))
        All_ratio(set_num) = 50;
    end
end  
OP_NUM = 8;    %������������
operator_name = {'CR', 'SV', 'MOP', 'SC', 'SCL', 'SM', 'SNS', 'WARP'};
operator_id = {'cr', 'sv', 'multiop', 'sc.', 'scl', 'sm', 'sns', 'warp'}; 

                                                                    
%%%%%%%%%%%
for set_num = 1:SET_NUM
    disp(['>>> #' num2str(set_num, '%03.0f') ' set --- ' PATH_NAME{set_num} '.']);
    % read the image set
    Path = [PATH_ROOT PATH_NAME{set_num} '\'];  %��������\��β
    file = dir([Path,'*.png']);    %����·���ļ����������ʽ���ļ����г�
     I2=im2double(imread([Path file(1).name]));    %ԭͼ
    retarget_name = zeros(OP_NUM,1);  %8��һ��
    
    %ѡ����������İ��ַ�����ͼƬ
    for i = 1:OP_NUM
        for j = 1 : size(file,1)      %����ԭͼ+��Сͼ������
            k1 = strfind(file(j).name, operator_id{i});  %ѡ���ض�ĳ��_id����
            if(All_ratio(set_num) == 75)
                k2 = strfind(file(j).name, '_0.75');
            elseif(All_ratio(set_num) == 50)
                k2 = strfind(file(j).name, '_0.50');
            end
            if( ~isempty(k1) && ~isempty(k2))   %~�ǣ�K1��K2��Ϊ�ǿյ�
                retarget_name(i) = j;    %�����������ķ����ź�
            end
        end
    end
   for op_num = 1:OP_NUM
        im_ret_set{op_num} =imread([Path file(retarget_name(op_num)).name]);%��Ӧ�ļ���  
   end
    %%%% �Ǹ�����׼
    
   for  op_num = 1:OP_NUM
       % I1=im2double(imread('D:\HYT\ABdata\BedRoom\BedRoom_0.75_cr.png'));
        %I2=im2double(imread('D:\HYT\Abdata\BedRoom\BedRoom.png'));
         I1 =im2double(im_ret_set{op_num});
        options.type='sd';  %FIR����Ӧ�˲���% Type of registration error used see registration_error.m
        options.centralgrad=false;   % Use fast forward instead of central error gradient
        Spacing=[4 4];  % b-spline grid spacing in x and y direction�����С��32*32
        [O_trans]=make_init_grid(Spacing,size(I1));  %�� % Make the Initial b-spline registration grid������ʼ��׼����,����Ǿ�������
        I1=double(I1); I2=double(I2); O_trans=double(O_trans);   % Convert all values to type double
         I2=imresize(I2,[size(I1)]);     %%%��ʱ 

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
        % fminsd�ҵ����������ĺ����ľֲ���Сֵ��ʹ����Ľ����½��Ż�  ��Դ˾��Ż� help fminsd��ϸ
        O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim); %��

        % Reshape O_trans from a vector to a matrix
        O_trans=reshape(O_trans,sizes);

        % Transform the input image with the found optimal grid.����bspline_transform
        Icor=bspline_transform(O_trans,I1,Spacing); %��

        % Make a (transformed) grid image
        Igrid=make_grid_image(Spacing,size(I1));  %��
        [Igrid,Tx,Ty,Tmap]=bspline_transform(O_trans,Igrid,Spacing);  %��jw[I,Tx,Ty,Tmap,Tz
        Igrid_all{set_num,op_num}=Igrid;
        Tx_all{set_num,op_num}= Tx;
        Ty_all{set_num,op_num}=Ty;
        Tmap_all{set_num,op_num}= Tmap;
%         figure,
%             imshow(Igrid); title('grid');
%             figure,imshow(Tmap);
    end
end



% Show the registration results
% figure,
% subplot(2,2,1), imshow(I1); title('input image 1');
% subplot(2,2,2), imshow(I2); title('input image 2');
% subplot(2,2,3), imshow(Icor); title('transformed image 1');
% subplot(2,2,4), imshow(Igrid); title('grid');
