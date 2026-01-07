function main()
%% Mesh settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.hmax = 0.03;%data.hmax = 0.03;
data.hmin = 0.005;%data.hmin = 0.005;
data.dif=100;
data.difc=5;
data.newmesh=2; % 1/2 new mesh/load mesh
%% Time settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t=0.034; % 25 fps = 0.04
t_0=0;
t_end=10;
theta = 0.55; % Scheme 0/0.5/1 Forward-Euler/Crank-Nicolson/Backward-Euler
data.post = 1; % 1/2 1-novideo/2-makevideo
%% Flow solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data.kyn=0.001; % kynematic viscosity
    data.rho=1; % hustota
    data.nu=data.kyn/data.rho; % dynamic
    data.vel=2; 
    % vel = 2 - SNS 78 ossen iteraci
    % vel = 2 - NNS 
%% Inlet function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SNS v(y) func spec 
 data.v_sns= @(y,t) 1.5*data.vel*4/0.1608.*(0.2+y).*(0.20-y);
 %% NNS v(y,t) func spec
 data.v_nns = @(y,t) 1.5*data.vel*4/0.1608.*(0.2+y).*(0.20-y).*...
    (0.5*(1-cos(pi/2*t))*(t<4)+(t>=8)*0.5*(1-cos(pi/2*t)));
 data.v_nns=data.v_sns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.problemtype = 2; %1/2 static SNS/nonstatic NNS
data.ns_tol = 1e-6; % ossen_precision
data.ns_max_it = 100;
data.print_ns_it = 2; % 1/2 1-no/2-yes
%% end settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data] = geometry_problem_DN(data);
%% Problem solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(data.problemtype==1)
    np=size(data.p,1);
        data.d_dir=compute_ddir(data);
    [uu,pp,it] = simulate_DirNeu_SNS(data);
    ux=uu(1:2:2*np-1); uy=uu(2:2:2*np);

    %% Plot actual results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Velocity %%%
    figure(41),clf,hold on
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),sqrt(ux.^2+uy.^2),'FaceColor','interp','LineStyle','none')
    colorbar
    view(2);
    axis equal
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);
    title('Velocity field');
    %% Pressure %%%
    figure(42),clf,hold on
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),pp,'FaceColor','interp','LineStyle','none')
    axis equal
    colorbar
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);
    view(2);
    title('Numeric pressure');
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(data.problemtype==2)
%% Time solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(data.post==2)
c=regexprep(regexprep(string(datetime),' ','_'),':','_');
v=VideoWriter('./Results/'+c+'.avi');
open(v);
end
%% Sets %%%
np=size(data.p,1);
ut10=zeros(np,1); ut20=ut10;
t_steps=t_0:delta_t:t_end;

for i=1:length(t_steps)
    fprintf('Frame %d/%d:\n',i,length(t_steps));
    ttt = t_steps(i);
    data.ttt=ttt;
    data.d_dir=compute_ddir(data);

    [uu,pp,it] = simulate_DirNeu_NNS(data,delta_t,theta,ut10,ut20);
    ux=uu(1:2:2*np-1); uy=uu(2:2:2*np); ut1=ux; ut2=uy;
    %% Plot actual results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Velocity %%%
    figure(41),clf,hold on
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),sqrt(ux.^2+uy.^2),'FaceColor','interp','LineStyle','none')
    colorbar
    view(2);
    axis equal
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);
    title('Velocity field');
    %% Pressure %%%
    figure(42),clf,hold on
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),pp,'FaceColor','interp','LineStyle','none')
    axis equal
    colorbar
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);
    view(2);
    title('Numeric pressure');
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(data.post==2)
    figure(40),clf,hold on
    subplot(2,1,1)    
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),sqrt(ux.^2+uy.^2),'FaceColor','interp','LineStyle','none')
    colorbar
    view(2);
    axis equal
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);  
    title(['Velocity field at time t = ',num2str(ttt, '%.3f.'), ' s']);
    %% Pressure %%%
    subplot(2,1,2)
    trisurf(data.t(:,1:3),data.p(:,1),data.p(:,2),pp,'FaceColor','interp','LineStyle','none')
    axis equal
    colorbar
    xlim([min(data.p(:,1)),max(data.p(:,1))]);
    ylim([min(data.p(:,2)),max(data.p(:,2))]);
    view(2);
    title(['Numeric pressure at the time t = ',num2str(ttt, '%.3f.'), ' s']);
       %% Uložení snímku %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        frame = getframe(gcf);
        for j=1:3 % size of video
         writeVideo(v,frame);
        end
    end
    %% Plot end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ut10=ut1; ut20=ut2;
end
if(data.post==2)
close(v);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Script done!");
end