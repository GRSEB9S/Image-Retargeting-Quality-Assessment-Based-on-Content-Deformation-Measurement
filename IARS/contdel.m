%过度集中删除度 
function delet_ratio=contdel(All_ratio,set_num,Func_aprox_Y)

[M,N]=size(Func_aprox_Y);
m=M-rem(double(M),16);
n=N-rem(double(N),16);
b_size=16;
step=8;

del_num=0;
for i=1:step:m-b_size
     for j=1:step:n-b_size
         black=numel(find(Func_aprox_Y(i:i+b_size,j:j+b_size)==-1));
         if(black/(b_size^2)>=0.9)
             del_num=del_num+1; 
         end
     end
end

delet_ratio=(del_num*b_size^2)*100/(M*N*All_ratio(set_num));
