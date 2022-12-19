% lsq_OCS
function res=lsq_BC(par)

bc_model=bc_optiloop(par);

%bcmod=[bc_model(end,3)+bc_model(end,10),bc_model(end,7)+bc_model(end,14)];
bcmod=bc_model(end,1:7)+bc_model(end,8:14);
totconc=sum(bcmod(:));
bcmodage=(bc_model(end,1:7).*bc_model(end,15:21)+bc_model(end,8:14).*bc_model(end,22:28))./totconc;

%bcopti=[2.2 2.0 2.2 2.6 2.0 1.2 1.2]; % max
%bcopti=[1.4 1.4 1.4 1.4 1.3 1.2 1.2]; %min
bcopti=[1.5 1.4 2.1 1.8 1.6 1.3 1.2];
bcoptiage=[8000 8000 4800 8000 15000 24000 24000];

diffconc=abs(bcopti-bcmod);
diffage=abs(bcmodage-bcoptiage);
%res=[diffconc/mean(diffconc), diffage/mean(diffage)];
res=diffconc;