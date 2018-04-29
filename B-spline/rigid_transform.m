function Iout=rigid_transform(Iin,M)
% Function rigid_transform, is a wrapper of the mex 
% rigid_transform_2d_double and rigid_transform_3d mex functions
%
% Iout = rigid_transform(Iin,M)
%
% inputs,
%   Iin :  Input image. 可以处理二维和三维的图
%   M : Transformation matrix
%
% outputs,
%   Iout: Output image
%
% Function is written by D.Kroon University of Twente (September 2008)
if(ndims(Iin)==2)     %输入图像是否为2维
    if(~isa(Iin,'double')), Iin=im2double(Iin); end    %判定是否为double型，不是的话转成double
	Iout=rigid_transform_2d_double(double(Iin),double(M));
else
    if(isa(Iin,'double'))
		Iout=rigid_transform_3d_double(double(Iin),double(M));
    else
		Iout=rigid_transform_3d_single(single(Iin),single(M));
    end
end
