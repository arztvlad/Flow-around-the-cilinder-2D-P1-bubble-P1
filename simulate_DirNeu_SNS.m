function [uu,pp,outer_it]=simulate_DirNeu_SNS(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: GENERATE PROBLEM (produces data structure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np=size(data.p,1);
w=zeros(2*np,1); % Leads to the Stokes problem in the first iteration
w(data.ind_dir)=data.d_dir;
nuu=size(data.ind_rest,1);

if(data.print_ns_it == 2)
fprintf('Oseen_iter |  err  \n')
fprintf('-----------------------------------------\n')
form='%4d | %8.2d \n';
end

reseni_ossen_0=zeros(nuu+np,1);
% Main loop
err=1; outer_it=0;
while err>data.ns_tol && outer_it<data.ns_max_it

    [A,Bu,Bl,E,b,c,~] = assemblyP1bP1DC_vec(data.p,data.t,data.nu,0,data.f,w);
    [A,Bu,Bl,E,b,c]=dirneu_BC(A,Bu,Bl,E,b,c,data.ind_dir,data.d_dir,data.p,[],[],[]);


    reseni_ossen=[A Bu'; Bl -E]\[b;c];
    err=norm(reseni_ossen-reseni_ossen_0)/(norm(reseni_ossen)+1);
    reseni_ossen_0=reseni_ossen;   
    w(data.ind_rest,1)=reseni_ossen(1:nuu);
  
    outer_it = outer_it+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(data.print_ns_it == 2)
    stat=[outer_it,err];
    fprintf(form,stat)
    end

end
uu=w;
pp=reseni_ossen(nuu+1:end);
fprintf('=========================================\n')
end