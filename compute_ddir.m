function [d_dir]=compute_ddir(data)

nd=size(data.ID,1);
d1=zeros(nd,1); d2=zeros(nd,1);
y=data.p(data.id.left,2);

if(data.problemtype==1)
    d1(1:length(data.id.left))=data.v_sns(y);
elseif(data.problemtype==2)
    d1(1:length(data.id.left))=data.v_nns(y,data.ttt);
end

d_dir(1:2:2*nd-1,1)=d1; d_dir(2:2:2*nd,1)=d2;
end