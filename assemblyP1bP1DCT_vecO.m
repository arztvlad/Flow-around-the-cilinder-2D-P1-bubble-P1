function [A,Bu,Bl,E,b,c,a]=assemblyP1bP1DCT_vecO(p,t,nu,f,w,delta_T,theta,u1,u2)

 np=size(p,1); alpha=1/delta_T;
 % 1. triangles area
 x21=p(t(:,2),1)-p(t(:,1),1); y12=p(t(:,1),2)-p(t(:,2),2);
 x32=p(t(:,3),1)-p(t(:,2),1); y23=p(t(:,2),2)-p(t(:,3),2);
 x13=p(t(:,1),1)-p(t(:,3),1); y31=p(t(:,3),2)-p(t(:,1),2);
 tarea=(x21.*y31-x13.*y12)/2;
 % 2. x^(T), y^(T), z^(T) and omega
 xt=[y23 y31 y12]; yt=[x32 x13 x21];
 
 nut=theta*nu./tarea;
 omega_12=(81/80)*nut.*(y23.*x32+y31.*x13+y12.*x21);
 omega_11=((81*alpha)/280)*tarea+(81/20)*nut.*(y23.^2-y31.*y12+0.5*(x32.^2-x13.*x21));
 omega_22=((81*alpha)/280)*tarea+(81/20)*nut.*(0.5*(y23.^2-y31.*y12)+x32.^2-x13.*x21);
 invomega=1./(omega_11.*omega_22-omega_12.^2);
 
am_12=omega_12.*invomega;
am_22=omega_22.*invomega;
am_11=omega_11.*invomega;

w1=w(1:2:end-1); w2=w(2:2:end);
w1t=(w1(t(:,1))+w1(t(:,2))+w1(t(:,3)))/3;
w2t=(w2(t(:,1))+w2(t(:,2))+w2(t(:,3)))/3;
cal=(theta/6)*[w1t.*xt(:,1)+w2t.*yt(:,1),w1t.*xt(:,2)+w2t.*yt(:,2),w1t.*xt(:,3)+w2t.*yt(:,3)];

u1t=(u1(t(:,1))+u1(t(:,2))+u1(t(:,3)))/3; 
u2t=(u2(t(:,1))+u2(t(:,2))+u2(t(:,3)))/3;
ttarea=2*tarea;
dxu1=((u1(t(:,2))-u1(t(:,1))).*y31+(u1(t(:,3))-u1(t(:,1))).*y12)./ttarea; 
dyu1=((u1(t(:,3))-u1(t(:,1))).*x21+(u1(t(:,2))-u1(t(:,1))).*x13)./ttarea;
dxu2=((u2(t(:,2))-u2(t(:,1))).*y31+(u2(t(:,3))-u2(t(:,1))).*y12)./ttarea; 
dyu2=((u2(t(:,3))-u2(t(:,1))).*x21+(u2(t(:,2))-u2(t(:,1))).*x13)./ttarea;

     bc1 = dxu1.*u1t+dyu1.*u2t;
     bc2 = dxu2.*u1t+dyu2.*u2t;

cl=(9/40)*(w1t.*xt+w2t.*yt);
zu =((3*alpha)/20)*tarea-theta*cl;
zl =((3*alpha)/20)*tarea+theta*cl;

f1=f(1:2:end-1);f2=f(2:2:end);


fh1=(f1(t(:,1))+f1(t(:,2))+f1(t(:,3)))/3; 
fh2=(f2(t(:,1))+f2(t(:,2))+f2(t(:,3)))/3; 

  f1t=tarea./3.*fh1+alpha*tarea./3.*u1t-(1-theta)*tarea./3.*bc1;
  f2t=tarea./3.*fh2+alpha*tarea./3.*u2t-(1-theta)*tarea./3.*bc2;  

  f1b=9/20*tarea.*fh1+alpha*9/20*tarea.*u1t-(1-theta)*9/20*tarea.*bc1;   
  f2b=9/20*tarea.*fh2+alpha*9/20*tarea.*u2t-(1-theta)*9/20*tarea.*bc2;   

% constants indemended on i
oyoi=(am_22.*f1b-am_12.*f2b);
oxoi=(am_11.*f2b-am_12.*f1b);
 
 A1=sparse(np,np); A2=sparse(np,np);
 A12=sparse(np,np);A21=sparse(np,np);
 E=sparse(np,np);
 Bu1=sparse(np,np); Bu2=sparse(np,np);
 Bl1=sparse(np,np); Bl2=sparse(np,np);
 b1=zeros(np,1); b2=zeros(np,1);
 c=zeros(np,1); 
a=zeros(np,1);
 for i=1:3
     a=a+sparse(t(:,i),1,tarea/3,np,1);
     % the right side vector components
     b1=b1+sparse(t(:,i),1,f1t -zu(:,i).*oyoi,np,1);
     b2=b2+sparse(t(:,i),1,f2t -zu(:,i).*oxoi,np,1);
     c=c-sparse(t(:,i),1,(9/40)*(xt(:,i).*oyoi+yt(:,i).*oxoi),np,1);
     % parts for matrices indemended on j
    zuom22=zu(:,i).*am_22;
    zuom11=zu(:,i).*am_11;
    zuom12=zu(:,i).*am_12;
    nuxt=(nut/4).*yt(:,i);
    nuyt=(nut/4).*xt(:,i);
    E1=(9/40)*(am_22.*xt(:,i)-am_12.*yt(:,i));
    E2=(9/40)*(am_11.*yt(:,i)-am_12.*xt(:,i));
     for j=1:3
         if(i==j)
             A1=A1+sparse(t(:,i),t(:,j),(alpha/6)*tarea+nuxt.*yt(:,j)+2*nuyt.*xt(:,j)+cal(:,j)-zuom22.*zl(:,j),np,np);
             A2=A2+sparse(t(:,i),t(:,j),(alpha/6)*tarea+2*nuxt.*yt(:,j)+nuyt.*xt(:,j)+cal(:,j)-zuom11.*zl(:,j),np,np);
         else
             A1=A1+sparse(t(:,i),t(:,j),(alpha/12)*tarea+nuxt.*yt(:,j)+2*nuyt.*xt(:,j)+cal(:,j)-zuom22.*zl(:,j),np,np);
             A2=A2+sparse(t(:,i),t(:,j),(alpha/12)*tarea+2*nuxt.*yt(:,j)+nuyt.*xt(:,j)+cal(:,j)-zuom11.*zl(:,j),np,np);
         end
         A12=A12+sparse(t(:,i),t(:,j),nuxt.*xt(:,j)+zuom12.*zl(:,j),np,np);
         A21=A21+sparse(t(:,i),t(:,j),nuyt.*yt(:,j)+zuom12.*zl(:,j),np,np);
         Bu1=Bu1-sparse(t(:,i),t(:,j),xt(:,j)/6+zu(:,j).*E1,np,np);
         Bu2=Bu2-sparse(t(:,i),t(:,j),yt(:,j)/6+zu(:,j).*E2,np,np);
         Bl1=Bl1-sparse(t(:,i),t(:,j),xt(:,j)/6+zl(:,j).*E1,np,np);
         Bl2=Bl2-sparse(t(:,i),t(:,j),yt(:,j)/6+zl(:,j).*E2,np,np);
         E=E+sparse(t(:,i),t(:,j),(9/40)*(E1.*xt(:,j)+E2.*yt(:,j)),np,np);
     end
 end
 t2(1:2:2*np-1)=1:np; t2(2:2:2*np)=np+(1:np);
 A=[A1 A12; A21 A2]; A=A(t2,t2);
 Bu=[Bu1 Bu2]; Bu=Bu(:,t2);
 Bl=[Bl1 Bl2]; Bl=Bl(:,t2);
 b=[b1;b2]; b=b(t2,1);

end