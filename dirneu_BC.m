function [A,Bu,Bl,E,b,c]=dirneu_BC(A,Bu,Bl,E,b,c,ind_dir,d,p,IN,s1,s2)

% Neumann boundary condition
if(~isempty(IN))
L=sqrt(sum((p(IN(:,1),:)-p(IN(:,2),:)).^2,2)); % length of the edges
for i=1:size(IN,1)
    eind=IN(i,1:2); % index of nodes on boundary edge
    b1=(s1(i,1)+s1(i,2))*L(i)/4;
    b2=(s2(i,1)+s2(i,2))*L(i)/4;
    v1ind=2*eind-1; v2ind=2*eind;
    b([v1ind,v2ind])=b([v1ind,v2ind])+[b1;b1;b2;b2];
end
end

% Dirichlet boundary condition
b=b-A(:,ind_dir)*d;
A(ind_dir,:)=[]; A(:,ind_dir)=[];
c=c-Bl(:,ind_dir)*d;
b(ind_dir)=[]; 
Bu(:,ind_dir)=[];
Bl(:,ind_dir)=[];


%d_p=zeros(length(p_dir),1);
%b=b-Bu(p_dir,:)'*d_p;


%c=c-Bl(:,p_dir)*d_p;
%  c(p_dir)=[]; 
% Bu(p_dir,:)=[];
%  Bl(p_dir,:)=[];
% E(p_dir,:)=[]; 
% E(:,p_dir)=[];

 % K = [A Bu'; Bl -E];
 % f = [b;c];
 % K*f;

 %c=c-Bl(:,p_dir)*d_p;

% c(p_dir)=10; 
% 
%  Bu(p_dir,:)=0;
%   Bl(p_dir,:)=0;
%  E(p_dir,:)=0; 
%  E(:,p_dir)=0;
% 
% E(sub2ind(size(E), p_dir, p_dir)) = 1;  % nastav diagonálu na 1 

end