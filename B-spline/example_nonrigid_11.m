% Example of nonrigid registration using steepest gradient optimizer使用最陡的梯度优化器的非刚性配准的示例
% Example is written by D.Kroon University of Twente (October 2008)
clear all; close all; 

%% path and other initial information
PATH_ROOT = 'D:\HYT\gray\'; % the path direct to the MIT dataset
load([PATH_ROOT 'subjData-ref_37.mat'])
subj_data = subjData.data;
SET_NUM = 37; % 37 set images
C1 = 1e-6;   
PATH_NAME = cell(SET_NUM,1);        %37行一列的空矩阵
All_ratio = zeros(SET_NUM,1);       %37行一列的零矩阵
for set_num = 1:SET_NUM            %循环37次    
    foo = subjData.datasetNames{set_num}; % 第{}个子图像的名字
    foo_loc = strfind(foo,'_0.');       %寻找foo中带有‘’中的位置
    PATH_NAME{set_num} =  foo(1:foo_loc-1);   %路径名字显示到'_0'前一位
    if(strfind(foo,'_0.75'))           %显示缩放比
        All_ratio(set_num) = 75;     
    elseif(strfind(foo,'_0.50'))
        All_ratio(set_num) = 50;
    end
end  
OP_NUM = 8;    %操作方法数量
operator_name = {'CR', 'SV', 'MOP', 'SC', 'SCL', 'SM', 'SNS', 'WARP'};
operator_id = {'cr', 'sv', 'multiop', 'sc.', 'scl', 'sm', 'sns', 'warp'}; 

                                                                    
%%%%%%%%%%%
for set_num = 1:SET_NUM
    disp(['>>> #' num2str(set_num, '%03.0f') ' set --- ' PATH_NAME{set_num} '.']);
    % read the image set
    Path = [PATH_ROOT PATH_NAME{set_num} '\'];  %连起来以\结尾
    file = dir([Path,'*.png']);    %将此路径文件夹中这个格式的文件名列出
     I2=im2double(imread([Path file(1).name]));    %原图
    retarget_name = zeros(OP_NUM,1);  %8行一列
    
    %选出测试所需的八种方法的图片
    for i = 1:OP_NUM
        for j = 1 : size(file,1)      %里面原图+缩小图的总数
            k1 = strfind(file(j).name, operator_id{i});  %选择特定某种_id方法
            if(All_ratio(set_num) == 75)
                k2 = strfind(file(j).name, '_0.75');
            elseif(All_ratio(set_num) == 50)
                k2 = strfind(file(j).name, '_0.50');
            end
            if( ~isempty(k1) && ~isempty(k2))   %~非，K1，K2均为非空的
                retarget_name(i) = j;    %给参与评估的方法排号
            end
        end
    end
   for op_num = 1:OP_NUM
        im_ret_set{op_num} =imread([Path file(retarget_name(op_num)).name]);%对应文件名  
   end
    %%%% 非刚性配准
    
   for  op_num = 1:OP_NUM
       % I1=im2double(imread('D:\HYT\ABdata\BedRoom\BedRoom_0.75_cr.png'));
        %I2=im2double(imread('D:\HYT\Abdata\BedRoom\BedRoom.png'));
         I1 =im2double(im_ret_set{op_num});
        options.type='sd';  %FIR自适应滤波器% Type of registration error used see registration_error.m
        options.centralgrad=false;   % Use fast forward instead of central error gradient
        Spacing=[4 4];  % b-spline grid spacing in x and y direction网格大小是32*32
        [O_trans]=make_init_grid(Spacing,size(I1));  %调 % Make the Initial b-spline registration grid创建初始配准网格,输出是均匀网格
        I1=double(I1); I2=double(I2); O_trans=double(O_trans);   % Convert all values to type double
         I2=imresize(I2,[size(I1)]);     %%%临时 

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
        % fminsd找到几个变量的函数的局部最小值，使用最陡的渐变下降优化  针对此句优化 help fminsd详细
        O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim); %调

        % Reshape O_trans from a vector to a matrix
        O_trans=reshape(O_trans,sizes);

        % Transform the input image with the found optimal grid.调用bspline_transform
        Icor=bspline_transform(O_trans,I1,Spacing); %调

        % Make a (transformed) grid image
        Igrid=make_grid_image(Spacing,size(I1));  %调
        [Igrid,Tx,Ty,Tmap]=bspline_transform(O_trans,Igrid,Spacing);  %调jw[I,Tx,Ty,Tmap,Tz
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
