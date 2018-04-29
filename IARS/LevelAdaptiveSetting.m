function [ levels ] = LevelAdaptiveSetting(img_org)
% Summary of this function goes here
%   Detailed explanation goes here
    [h_org, w_org] = size(img_org);
    if(h_org < w_org)
        levels =  ceil(log2(w_org/10));  %大于等于一个近似的整数
    else
        levels =  ceil(log2(h_org/10));
    end

end
%levels={log[max(w,h)/10]}
