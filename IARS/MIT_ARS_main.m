% The code is for research purpose only.
% experimental test on MIT RetargetMe dataset, 
% Code by Yabin Zhang, for the ICASSP paper
% [1] Yabin Zhang, Weisi Lin, Xinfeng Zhang, Yuming Fang, Leida Li.
% Aspect Ratio Similarity (ARS) for Image Retargeting Quality Assessment. ICASSP 2016
% [2] Yabin Zhang, Yuming Fang, Weisi Lin, Xinfeng Zhang, Leida Li.
% Backward Registration-Based Aspect Ratio Similarity for Image Retargeting Quality Assessment.
% Transaction on Image Processing (TIP), 2016.

% To run the code
% prepare the RetargetMe dataset and correct the root path

%clear all; 
%% path and other initial information
PATH_ROOT = 'E:\HYT\���ݿ�����\'; % the path direct to the MIT dataset
% load the subjective data (put the subjData at the same path)
load([PATH_ROOT 'subjData-ref_37.mat'])
subj_data = subjData.data;
SET_NUM = 1;   % 37 set images
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

%% Backward Registration

disp(' -------------------------------------------------------------------------------'); %''����ʾ����
disp('>>> start backward registration ...');
for set_num = 1:SET_NUM
    disp(['>>> #' num2str(set_num, '%03.0f') ' set --- ' PATH_NAME{set_num} '.']);
    % read the image set
    path = [PATH_ROOT PATH_NAME{set_num} '\'];  %��������\��β
    file = dir([path,'*.png']);    %����·���ļ����������ʽ���ļ����г�
    im_org =  imread([path file(1).name]);
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
        im_ret_set{op_num} = imread([path file(retarget_name(op_num)).name]);%��Ӧ�ļ���
    end
%  smap =  imread(['.\MIT_smap\' PATH_NAME{set_num} '_smap.png']); %��Ӧԭͼ����ͼ
     smap =  imread(['..\MIT_smap\' PATH_NAME{set_num} '_smap.png']);
    All_img_org{set_num} = im_org  ;        %һ��37��
    All_smap{set_num} = smap;                %������ͼ��š�
    for op_num = 1:OP_NUM
        disp(['  ---+ #' num2str(op_num, '%02.0f') ' retargeted image: ' operator_name{op_num}]);        
        im_ret = im_ret_set{op_num};   %��ǰ��ŵ�ͼ
        %[foo_XX, foo_YY] = BWRegistration(im_org, im_ret);     %���ú�����׼
        All_img_ret{set_num,op_num} = im_ret;
        %All_XX{set_num,op_num} = foo_XX;
        %All_YY{set_num,op_num} = foo_YY;
    end
end
   
%% ARS evaluation
BLK_SIZE = 16;%16
ALPHA = 0.30;

disp(' -------------------------------------------------------------------------------');
% reveal the forward retargeting information��ʾ
All_BLK_changes = cell(SET_NUM,OP_NUM);  %37*8
All_BLK_sal = cell(SET_NUM,1);
for set_num = 1:SET_NUM
    im_org = All_img_org{set_num};
    smap = double(All_smap{set_num});
    disp(['->-> #' num2str(set_num, '%02.0f') ' [' PATH_NAME{set_num} ']  image set evaluating ...']);
    [height_org, width_org,~] = size(im_org);
    blk_h = floor(height_org/BLK_SIZE); blk_w = floor(width_org/BLK_SIZE);     %floor �ͽ�����
    blk_sal_org = zeros(blk_h, blk_w);
    smap = smap/sum(smap(:));  %��һ��
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
        All_Func_Y{set_num,op_num}=Func_aprox_Y;   %����ɾ��
        delet_ratio=contdel(All_ratio,set_num,Func_aprox_Y);
        All_delet_ratio(set_num,op_num)=delet_ratio;
        [Block_change_info,dist_ratio] = ReTransBLK(im_org, Func_aprox_X, Func_aprox_Y, BLK_SIZE);
        All_BLK_changes{set_num, op_num} = Block_change_info;
         All_BLK_ratio{set_num, op_num}=dist_ratio;%
    end
end

score_ARS = zeros(SET_NUM,OP_NUM);
for set_num = 1:SET_NUM
    % collect the data for the set
    CBlock_sal_org = All_BLK_sal{set_num};
    [blk_h, blk_w] = size(CBlock_sal_org);
    ARS = zeros(blk_h, blk_w); 
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
                 %m_ratio=CBlock_ratio(bi,bj);%
                 %ARS(bi, bj) = exp( -ALPHA*(m_ratio-1).^2)*((2*w_ratio*h_ratio+C1)/(w_ratio^2+h_ratio^2+C1))^0.5;
                ARS(bi, bj) = exp( -ALPHA*(m_ratio*10/All_ratio(set_num)-1).^2)*((2*w_ratio*h_ratio+C1)/(w_ratio^2+h_ratio^2+C1))^0.16;
            end
        end  
        foo_score = CBlock_sal_org.*ARS;
        score_ARS(set_num, op_num) =sum(foo_score(:));
    end    
end

% KRCC evaluation
KRCC = zeros(SET_NUM,1);
for set_num = 1:SET_NUM
    KRCC(set_num) = getKLCorr(subj_data(set_num,:), score_ARS(set_num,:));  %�ϵ¶�����ʽ
end

disp(' -------------------------------------------------------------------------------');
disp(['Block size -- '  num2str(BLK_SIZE) '; coefficient alpha = ' num2str(ALPHA)]);
disp(['mean KRCC = ' num2str(mean(KRCC), '%0.3f')]);
disp(['std KRCC = ' num2str(std(KRCC), '%0.3f')]);


