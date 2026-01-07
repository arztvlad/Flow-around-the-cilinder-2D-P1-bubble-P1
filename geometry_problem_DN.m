function [data]=geometry_problem_DN(data)

  addpath('MESH2D')

  if(data.newmesh==1)
[p,etri,t,tnum]=rect_hol_mesh(data.hmax,data.hmin,data.dif,data.difc);
save('mesh_SNS','p','etri','t','tnum');
  end
load('mesh_SNS','p','etri','t');
data.p=p;data.t=t; np=size(data.p,1);
idetri=unique(etri);
data.idetri=idetri;
[id.left, id.right, id.upper, id.below, id.circ,id.circ_mask] = make_edge(data.p,idetri);
[in.left, in.right, in.upper, in.below, in.circ] = make_edges(data.p,id);
data.id=id;
data.in=in;
data.p_u_idx = find(abs(p(:,2) - 0.2) < 1e-12);
data.p_b_idx = find(abs(p(:,2) + 0.2) < 1e-12);

% 2.1) Boundary nodes - for Dirichlet BC
% indices in vectors
ID= [id.left;id.circ;id.upper;id.below];
nd=size(ID,1);
data.ind_dir=reshape([2*ID'-1;2*ID'],2*nd,1);
data.ind_rest=setdiff([1:2*np]',data.ind_dir);
data.ID=ID;
%% 2.2) Boundary nodes - for Neumann BC [p1 p2]
data.IN=in.right;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SECTION 2: PDE PROBLEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Right hand-side function
% Experimental solution (analytic for Dirichlet and Neumann BC on [0,1]x[0,1])
data.f(1:2:2*np-1,1)=zeros(np,1); data.f(2:2:2*np,1)=zeros(np,1);
%% 3.1) Boundary conditions - Dirichlet
% d1=zeros(nd,1); d2=zeros(nd,1);
% y=data.p(data.id.left,2);
% d1(1:length(data.id.left))=(6/0.1608)*data.vel.*(0.2+y).*(0.20-y);
% d1(1:length(data.id.left))=data.v_sns(y);
% 
% d_dir(1:2:2*nd-1,1)=d1; d_dir(2:2:2*nd,1)=d2;
% data.d_dir=d_dir;

%% 3.2) Boundary conditions - Neumann
% s1=zeros(size(in.right,1),1); s2=s1;
% data.s1=[s1 s1]; data.s2=[s2 s2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vert,etri,tria,tnum]=rect_hol_mesh(hmax,hmin,dif,difc)
%addpath('mesh2D');
libpath();
%DEMO8 explore impact of "hill-climbing" mesh optimisations.

%---------------------------------------------- create geom.
    node = [0,-0.2;
        2.2,-0.2;
        2.2,0.2;
        0,0.2];
    edge = [1,2;
        2,3;
        3,4;
        4,1];

    adel = 2.*pi/64;%64
    amin = 0;
    amax = 2.*pi-adel;

    xcir = 0.2+0.05*cos(amin:adel:amax)';
    ycir = 0.0+0.05*sin(amin:adel:amax)';
    ncir = [xcir,ycir] ;
    numc = size(ncir,1);

    ecir(:,1) = [(1:numc-1)';numc];
    ecir(:,2) = [(2:numc-0)';1];

    ecir = ecir+size(node,1);
    edge = [edge; ecir];
    node = [node; ncir];

    hfun = @hfun8b;
%---------------------------------------------- do mesh-gen.

   
   [vert,etri,tria,tnum] = refine2(node,edge,[],[],hfun);
    
%---------------------------------------------- do mesh-opt.
 [vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum) ;
    
    figure(21),clf,hold on;
    subplot(2,1,1);
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
        num2str(size(tria,1))]) ;
        
        subplot(2,1,2);
    patch('faces',tria(:,1:3),'vertices',vert , ...
        'facevertexcdata' , hfun(vert), ...
        'facecolor','interp', ...
        'edgecolor','none') ;
    hold on; axis image off;
    title('MESH-SIZE function.');
   
    hvrt = feval(hfun,vert) ;
    
    tricost(vert,etri,tria,tnum,hvrt) ;
           
    drawnow;
    % 
     set(figure(21),'units','normalized');%, ...
      %   'position',[.05,.50,.30,.35]) ;
    % set(figure(2),'units','normalized', ...
    %     'position',[.35,.50,.30,.35]) ; 
     set(figure(22),'units','normalized');%, ...
       %  'position',[.05,.05,.30,.35]) ;  
    %rmpath('mesh2D');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hfun] = hfun8b(test)
       hup=exp(-dif*(test(:,2)+0.2).^2);
    hbot=exp(-dif*(test(:,2)-0.2).^2);
    hleft=exp(-dif*(test(:,1)-0).^2);
    hright=exp(-dif*(test(:,1)-2.2).^2);

    hcir = exp(-difc*(test(:,1)-0.2).^2-difc*(test(:,2)-0.0).^2 ) ;
    hedg=max(max(max(max(hup,hbot),hleft),hright),hcir);
    hfun = hmax - (hmax-hmin) * hedg  ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
        rmpath('MESH2D')
end