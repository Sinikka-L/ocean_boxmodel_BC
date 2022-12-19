%% BC model MAIN
function [bc]=bc_optiloop(par)
addpath('./external_functions')
addpath('./cbrewer')
%% set up model scenario
% Earth parameters (number of domains, pools, etc.)
PE = setup_model;
% Scenario parameters (CO2 emissions, model coupling, etc.)
PS = setup_run(PE);
% Land parameters (carbon turnover times, etc.)
%PL = setup_land(PE,PS);
% Ocean parameters (circulation and biological rates, temperatures, etc.)
PO = setup_ocean(PE,PS);

%% set an initial state

% package initial values into matrix  
%   This makes it easier to identify domains (rows) and pools (columns)
y0mat = PE.m0; % start with dummy matrix
y0mat(PE.Ibc1)=1; % put temperature values into matrix
y0mat(PE.Ibc2)=1; % ditto ocean nutrient, etc.
y0mat(PE.Iage1)=0;
y0mat(PE.Iage2)=0;
y0mat(PE.Iwage)=0;
%y0mat(PE.Iatm)=Cat_o;

%% Integrate forward

% convert matrix of initial values to a vector (required for ode solver)
y0=y0mat(PE.Ires); 

% set time for integration
trun=[0:PE.spery:PS.ytot*PE.spery];

% use matlab's ode solver to run model forward in time
%   first argument: pointer to carbon_climate_derivs.m, the guts of the model
%   second argument: start and end times [in seconds]
%   third argument: initial condition for all pools
options = odeset('abstol', 1e-18, 'reltol', 1e-13,...
    'NonNegative', ones(size(y0)));

% make structure for all model parametrs in PD
% convert all to day!!!!
%PD.bc1sink = [0.16 0.25 0.34 0.26 0 0 0]*1/800*120; %per year (but light intensity needs to be in there?)
%PD.bc1sink = [0.16 0.25 0.34 0.26 0 0 0]*1/800/365; % per day (but light intensity needs to be in there?)
%PD.bc1sink = [0.035	0.042	0.862	0.0605 0 0 0]*70/800/365; % 164=mean W/m2 radiation
%PD.bc1sink = [0.035	0.042	0.862	0.0605 0 0 0]*3./PO.V*10^12/12*1000/365; % 3 Tg/yr Stubbins
%PD.bc1sink = [0.035	0.042	0.862	0.0605 0 0 0]*9./PO.V*10^12/12*1000/365; % 3 Tg/yr Stubbins
% halt stop das ist eine first order sink!!!!!! nix mit 3 Tg verteilen!!!
PD.bc1sink = [0.035	0.042	0.862	0.0605 0 0 0]*par(1)/800/365 ... % UV s
           +[1 1 1 1 1 1 1]*1/par(3)/365;% 3 Tg/yr Stubbins, 10 and 3000 worked

PD.bc2sink = [0.035	0.042	0.862	0.0605 0 0 0]*par(2)/800/365 ... % UV s
           +[1 1 1 1 1 1 1]*1/par(4)/365;% 
       
PD.Svec=zeros(PE.nb,2);
partfac=par(5);
BCinputriver=[0 0 16.2 1.8 0 0 0]*partfac; % Tg/yr
BCinputsediment=[0 0 0 0 0 0 0]; % Tg/yr 10
BCinputaerosols=[0 0 0 0 0 0 0]; % Tg/yr
BCinput=BCinputriver+BCinputsediment+BCinputaerosols;
BCinputriver2=[0 0 16.2 1.8 0 0 0]*(1-partfac); % Tg/yr
BCinputsediment2=[0 0 0 0 0 0 par(6)]; % Tg/yr 10
BCinputaerosols2=[0 0 0 0 0 0 0]; % Tg/yr
BCinput2=BCinputriver2+BCinputsediment2+BCinputaerosols2;
%frac=0.8;
%PD.Svec(:,1)=frac*BCinput./PO.V*10^12/12*1000/365; %mmol/d/m3 % river input
%PD.Svec(:,2)=(1-frac)*BCinput./PO.V*10^12/12*1000/365;
PD.Svec(:,1)=BCinput./PO.V*10^12/12*1000/365; %mmol/d/m3 % river input
PD.Svec(:,2)=BCinput2./PO.V*10^12/12*1000/365;

ageinputriver1=[0 0 0 0 0 0 0]*365;
ageinputsediment1=[0 0 0 0 0 0 0]*365;
ageinputaerosol1=[0 0 0 0 0 0 0]*365;
BCagein1=ageinputriver1+ageinputsediment1+ageinputaerosol1;
ageinputriver2=[0 0 0 0 0 0 0]*365;
ageinputsediment2=[0 0 0 0 0 0 24000]*365; % 24000
ageinputaerosol2=[0 0 0 0 0 0 0]*365;
BCagein2=ageinputriver2+ageinputsediment2+ageinputaerosol2;
PD.Svecage(:,1)=BCagein1'; % mmol/C: age*input rate
PD.Svecage(:,2)=BCagein2'; % mmol/C: age*input rate
PD.wagedt=[0 0 0 0 1 1 1];


[t,y]=ode15s(@carbon_BC_derivs,trun,y0', options, PE, PO, PD); 

bc=y;
end