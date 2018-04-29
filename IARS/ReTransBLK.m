function [Block_change_info,dist_ratio] = ReTransBLK(img_org, Func_aprox_X, Func_aprox_Y, BLOCK_SIZE)
% Summary of this function goes here
%   Detailed explanation goes here
    
    [height_org, width_org,~] = size(img_org);
   % [height_org, width_org] = size(img_org);
    block_h = floor(height_org/BLOCK_SIZE);
    block_w = floor(width_org/BLOCK_SIZE);

    Block_change_info = zeros(block_h, block_w, 7);
    for bi = 1:block_h
        for bj = 1:block_w           
            top_h = (bi-1)*BLOCK_SIZE+1; top_w = (bj-1)*BLOCK_SIZE+1;               
            CBlock_Func_aprox_X = Func_aprox_X(top_h:(top_h+BLOCK_SIZE-1), ...
                top_w:(top_w+BLOCK_SIZE-1));
            CBlock_Func_aprox_Y = Func_aprox_Y(top_h:(top_h+BLOCK_SIZE-1), ...
                top_w:(top_w+BLOCK_SIZE-1));
            
            CSet_Ret = [CBlock_Func_aprox_X(:)...
                CBlock_Func_aprox_Y(:)]; 
             dist_ratio(bi, bj)= numel(find(CSet_Ret(:,1)==-1))/(BLOCK_SIZE*BLOCK_SIZE);
            % the point set in ret img mapped from the regular block in org img
            CSet_Ret(CSet_Ret(:,1) == -1,:) = [];
            CSet_Ret(CSet_Ret(:,2) == -1,:) = [];
            
            if(isempty(CSet_Ret)) % measure the point set
                X_MAX_ret = 0; X_MIN_ret = 0; 
                Y_MAX_ret = 0; Y_MIN_ret = 0;
            else
                temp1=CSet_Ret(:,1); %
                temp2=CSet_Ret(:,2);%
                if  (size(temp1,1)>4&size(temp2,1)>4)%
                
                        [~,ind1]=sort(temp1);%
                        X_MAX_ret =temp1(ind1(end-3)); X_MIN_ret = min(CSet_Ret(:,1)); %

                        [~,ind2]=sort(temp2);%
                        Y_MAX_ret =temp2(ind2(end-3)); Y_MIN_ret = min(CSet_Ret(:,2));%
                else%
                        X_MAX_ret = max(CSet_Ret(:,1)); X_MIN_ret = min(CSet_Ret(:,1)); 
                        Y_MAX_ret = max(CSet_Ret(:,2)); Y_MIN_ret = min(CSet_Ret(:,2));
                end%
            end
            Block_change_info(bi, bj,1) = (X_MAX_ret-X_MIN_ret+1);
            Block_change_info(bi, bj,2) = (Y_MAX_ret-Y_MIN_ret+1);
        end
    end  
    
end

