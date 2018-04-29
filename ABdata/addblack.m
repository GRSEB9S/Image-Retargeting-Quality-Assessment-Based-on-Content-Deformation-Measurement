
 load ('E:\HYT\NRID\Yuming_Sal\data\All_XX_n.mat');
 load ('E:\HYT\NRID\Yuming_Sal\data\All_YY_n.mat');
 load ('E:\HYT\NRID\NRIDdata\All_img_org_n.mat');
 load ('E:\HYT\NRID\NRIDdata\All_img_ret_n.mat');
 load ('E:\HYT\NRID\NRIDdata\subjData-ref_35.mat');
SET_NUM = 35;
OP_NUM = 5;    %methods number
operator_name = { 'MO', 'SCL', 'SC', 'SM', 'WARP'};
operator_id ={'multiop', 'scal', 'seam', 'shif', 'warp'};   

img_dir='E:\HYT\NRID\ABdata_N\';

for set_num = 1:SET_NUM
    img_name=subjData.datasetNames{set_num};
    im_org =rgb2gray(All_img_org{set_num});
    [h,w]=size(im_org);
    disp(['->-> #' num2str(set_num, '%02.0f') ]);
  
    for op_num = 1:OP_NUM      
        im_ret = rgb2gray(All_img_ret{set_num, op_num});
        [H,W]=size(im_ret);
        XX = All_XX{set_num, op_num}; YY = All_YY{set_num, op_num};
         im_ad=im_org;
        z=zeros(h,w);  %black
        for i_r=1:H
           for j_r=1:W
              a=All_YY{set_num, op_num}(i_r,j_r);
              b=All_XX{set_num, op_num}(i_r,j_r);
              z(a,b)=1;        %white      
           end
        end
        for  a=1:h
            for b=1:w
                if     z(a,b)==0;
                     im_ad(a,b)=1;
                end
            end
        end
       imwrite(im_ad,[img_dir img_name '_'  operator_id{op_num} '.jpg']); 
    end
end
            
    
 