%%
% the demo is the main function
% To run the code
% prepare the RetargetMe dataset and correct the root path
% clear all; 
%% path and other initial information
PATH_ROOT = 'E:\HYT\数据库已用\'; % the path direct to the MIT dataset
load([PATH_ROOT 'subjData-ref_37.mat']) % load the subjective data (put the subjData at the same path)
subj_data = subjData.data;
SET_NUM = 20;   % 37 set images
C1 = 1e-6;   
PATH_NAME = cell(SET_NUM,1);       
All_ratio = zeros(SET_NUM,1);       
for set_num = 1:SET_NUM          
    foo = subjData.datasetNames{set_num}; 
    foo_loc = strfind(foo,'_0.');       
    PATH_NAME{set_num} =  foo(1:foo_loc-1);   
    if(strfind(foo,'_0.75'))          
        All_ratio(set_num) = 75;     
    elseif(strfind(foo,'_0.50'))
        All_ratio(set_num) = 50;
    end
end  
OP_NUM = 8;    
operator_name = {'CR', 'SV', 'MOP', 'SC', 'SCL', 'SM', 'SNS', 'WARP'};
operator_id = {'cr', 'sv', 'multiop', 'sc.', 'scl', 'sm', 'sns', 'warp'};    

%% Backward Registration
disp(' -------------------------------------------------------------------------------'); 
for set_num =20:SET_NUM
    disp(['>>> #' num2str(set_num, '%03.0f') ' set --- ' PATH_NAME{set_num} '.']);
    % read the image set
    path = [PATH_ROOT PATH_NAME{set_num} '\'];  
    path2=['E:\HYT\代码整理小论文\ABdata\'  PATH_NAME{set_num} '\'];
    file = dir([path,'*.png']);   
    im_org =  imread([path file(1).name]);
    I2=rgb2gray(im_org);
    retarget_name = zeros(OP_NUM,1);  
    
    for i = 1:OP_NUM
        for j = 1 : size(file,1)      
            k1 = strfind(file(j).name, operator_id{i});  
            if(All_ratio(set_num) == 75)
                k2 = strfind(file(j).name, '_0.75');
            elseif(All_ratio(set_num) == 50)
                k2 = strfind(file(j).name, '_0.50');
            end
            if( ~isempty(k1) && ~isempty(k2))  
                retarget_name(i) = j;    
            end
        end
    end
    for op_num = 1:OP_NUM
        im_ret_set{op_num} = imread([path file(retarget_name(op_num)).name]);
        im_ret_rg{op_num} = imread([path2 file(retarget_name(op_num)).name]);%后向配准预测位置图
    end
    smap =  imread(['E:\HYT\代码整理小论文\IARS\MIT_smap\' PATH_NAME{set_num} '_smap.png']);
    All_img_org{set_num} = im_org  ;      
    All_smap{set_num} = smap;                
    for op_num = 1:OP_NUM
        disp(['  ---+ #' num2str(op_num, '%02.0f') ' retargeted image: ' operator_name{op_num}]);  
        im_ret = im_ret_set{op_num};   
        [foo_XX, foo_YY] = BWRegistration(im_org, im_ret);     %调用后向配准
        All_img_ret{set_num,op_num} = im_ret;
        All_XX{set_num,op_num} = foo_XX;
        All_YY{set_num,op_num} = foo_YY;
    end
end
   
%% evaluation
BLK_SIZE = 16;
ALPHA = 0.30;
disp(' -------------------------------------------------------------------------------');
% reveal the forward retargeting information
All_BLK_changes = cell(SET_NUM,OP_NUM);  
All_BLK_sal = cell(SET_NUM,1);
for set_num = 20:SET_NUM
    im_org = All_img_org{set_num};
    smap = double(All_smap{set_num});
    disp(['->-> #' num2str(set_num, '%02.0f') ' [' PATH_NAME{set_num} ']  image set evaluating ...']);
    [height_org, width_org,~] = size(im_org);
    blk_h = floor(height_org/BLK_SIZE); blk_w = floor(width_org/BLK_SIZE);    
    blk_sal_org = zeros(blk_h, blk_w);
    smap = smap/sum(smap(:));  %归一化
    for bi = 1:blk_h
        for bj = 1:blk_w
            top_h = (bi-1)*BLK_SIZE+1; top_w = (bj-1)*BLK_SIZE+1;
            CBlock_sal = smap(top_h:(top_h+BLK_SIZE-1), ...
                top_w:(top_w+BLK_SIZE-1));
            blk_sal_org(bi, bj) = sum(sum(CBlock_sal));
        end
    end
    All_BLK_sal{set_num} = blk_sal_org;
    for op_num = 1:OP_NUM      
        im_ret = All_img_ret{set_num, op_num};
        XX = All_XX{set_num, op_num}; YY = All_YY{set_num, op_num};
        [Func_aprox_X,   Func_aprox_Y] = ReforumlatedMapping(im_org, XX, YY);
        All_Func_Y{set_num,op_num}=Func_aprox_Y; 
        
        %Concentrated deletion ratio
        delet_ratio=contdel(All_ratio,set_num,Func_aprox_Y);
        All_delet_ratio(set_num,op_num)=delet_ratio;
        [Block_change_info,dist_ratio] = ReTransBLK(im_org, Func_aprox_X, Func_aprox_Y, BLK_SIZE);
        All_BLK_changes{set_num, op_num} = Block_change_info;
         All_BLK_ratio{set_num, op_num}=dist_ratio;
         
         %%foreground retention ratio
         img_rbd{set_num} = imread(['E:\HYT\MIT_rbd\rbd\' PATH_NAME{set_num} '_rgb.png']); %
         img_ab=cell2mat(im_ret_rg(op_num));   
         rat=forebi(img_rbd{set_num},img_ab,All_ratio(set_num));
         ratio_bi(set_num,op_num)=rat;
         
         % B-spline elastic registration
        I1 =im2double(im_ret_rg{op_num});
        options.type='sd'; % Type of registration error used see registration_error.m
        options.centralgrad=false;   % Use fast forward instead of central error gradient
        Spacing=[16 16];  
        [O_trans]=make_init_grid(Spacing,size(I1));  % Make the Initial b-spline registration grid
        I1=double(I1); I2=double(I2); O_trans=double(O_trans);   % Convert all values to type double
        I2=imresize(I2,[size(I1)]);     %%

        I1s=imfilter(I1,fspecial('gaussian')); % Smooth both images for faster registration
        I2s=imfilter(I2,fspecial('gaussian'));
     
        optim=struct('Display','iter','GradObj','on','MaxIter',20,'DiffMinChange',0.1,'DiffMaxChange',1);   % Optimizer parameters      
        sizes=size(O_trans); O_trans=O_trans(:);  % Reshape O_trans from a matrix to a vector.
        O_trans = fminsd(@(x)bspline_registration_gradient(x,sizes,Spacing,I1s,I2s,options),O_trans,optim); %Start the b-spline nonrigid registration optimizer  
        O_trans=reshape(O_trans,sizes);% Reshape O_trans from a vector to a matrix 
        Icor=bspline_transform(O_trans,I1,Spacing); %调% Transform the input image with the found optimal grid

        % Make a (transformed) grid image
        Igrid=make_grid_image(Spacing,size(I1));  %调
        [Igrid,Tx,Ty,Tmap]=bspline_transform(O_trans,Igrid,Spacing);  
        Igrid_all{set_num,op_num}=Igrid;
        Tx_all{set_num,op_num}= Tx;
        Ty_all{set_num,op_num}=Ty;
        Tmap_all{set_num,op_num}= Tmap;
        
         figure,imshow(Igrid); title('grid');
         Varsum_b=varsal(SET_NUM,OP_NUM,Tx,Ty,img_rbd);
              
    end
end

%%
%IARS evaluation
score_IARS = zeros(SET_NUM,OP_NUM);

    % collect the data for the set
    CBlock_sal_org = All_BLK_sal{set_num};
    [blk_h, blk_w] = size(CBlock_sal_org);
    IARS = zeros(blk_h, blk_w); 
    for op_num = 1:OP_NUM
        CBlock_change_info = All_BLK_changes{set_num, op_num};
        CBlock_ratio=All_BLK_ratio{set_num, op_num};%
        CBlock_info_w = CBlock_change_info(:,:,1);
        CBlock_info_h = CBlock_change_info(:,:,2);
        % compute the distortion for each op_num
        for bi = 1:blk_h
            for bj = 1:blk_w
                w_ratio = ( CBlock_info_w(bi, bj) )/BLK_SIZE;
                h_ratio = ( CBlock_info_h(bi, bj) )/BLK_SIZE;
                m_ratio = (w_ratio + h_ratio)/2;
                IARS(bi, bj) = exp( -ALPHA*(m_ratio*10/All_ratio(set_num)-1).^2)*((2*w_ratio*h_ratio+C1)/(w_ratio^2+h_ratio^2+C1))^0.16;
            end
        end  
        foo_score = CBlock_sal_org.*IARS;
        score_IARS(set_num, op_num) =sum(foo_score(:));
    end    

%%
final_score= score_IARS-(Varsum_b+All_BLK_ratio*0.01)*0.04+ratio_bi*0.02 ;

% KRCC evaluation
KRCC= getKLCorr(subj_data(set_num,:),final_score);
disp(['the score of car1 -- KRCC='  KRCC]);
