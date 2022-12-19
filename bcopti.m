%%% optimization
clc
%clear all
close all

options=optimset('lsqnonlin');
options.TolX=1e-11;
options.TolFun=1e-11;
options.Display='iter';
options.Jacobian='off';
options.Algorithm='levenberg-marquardt';
lb=[10,10,50,10000,0.1,0.1];
ub=[100,100,1000,5000000,1,3];


%par=[5*10-7, 4*10^-7, 2*10^-7, 1.5*10^-7, 9*10^-8, 2*10^-8]; % p1-5: different wavelengths apparent quantum yield
par=[10,10,100,100000,0.7,1.5];
[phat,resnorm,residual,exitflag,output,lambda,Jacob]=lsqnonlin(@lsq_BC,par, lb, ub, options);
ci=nlparci(phat, residual, 'jacobian', Jacob);  

%%
y=bc_optiloop(phat);

%% ziel range stock
%surfBC=sum(y(end,1:4).*PO.V(1:4))*12/10^15/1000 % Pg
%subsurfBC=sum(y(end,5:7).*PO.V(5:7))*12/10^15/1000 % Pg
surfBCmin=min(y(end,1:4)+y(end,8:11)) % 然 = mmol/m3
surfBCmax=max(y(end,1:4)+y(end,8:11)) % 然 = mmol/m3
subsurfBCmin=min(y(end,5:7)+y(end,12:14)) % 然 = mmol/m3
subsurfBCmax=max(y(end,5:7)+y(end,12:14)) % 然 = mmol/m3
agemin=min(mean([y(end,15:21);y(end,22:28)]))/365
agemax=max(mean([y(end,15:21);y(end,22:28)]))/365

% plot BC concentration
plot_boxes(y(end,1:7)+y(end,8:14),'clabel','BC conc. [mmol m^{-3}]');
% plot BC age
sumbctot=y(end,1:7)+y(end,8:14);
agemeanweighted=[y(end,15:21).*y(end,1:7)./sumbctot;y(end,22:28).*y(end,8:14)./sumbctot];
plot_boxes(sum(agemeanweighted)/365,'clabel','BC age [years]');
