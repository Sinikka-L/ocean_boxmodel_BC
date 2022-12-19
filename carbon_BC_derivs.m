function [dydt] = carbon_BC_derivs(t,y, PE, PO, PD)
%% Create local copies of state variables

BC1loc = y(PE.Jbc1)'; % Black Carbon 1 [mmol/m³]
BC2loc = y(PE.Jbc2)'; % Black Carbon 2 [mmol/m³]
BC1ageloc = y(PE.Jage1)'; % Radiocarbon age BC 1 [d]
BC2ageloc = y(PE.Jage2)'; % Radiocarbon age BC 2[d]
wageloc = y(PE.Jwage)'; % water masse age [d]


%% update basic quantities

Psi   = PO.Psi_o; % [1/d]


%% DOM ecosystem model - no temperature dependance [basic unit: mmolC/m³/d]


% compute first order sinks
BC1 = -BC1loc.*PD.bc1sink; % source Svec is added later
BC2 = -BC2loc.*PD.bc2sink;



%% Compute Tendencies - should have units mmol/d /m3

% tendency = physics + (sink, negative) + DBC-River+Sediment + DBC particular
% dissolved river + aerosols
dBC1dt = (Psi*BC1loc' + BC1') + PD.Svec(:,1)+PD.Svecriver(:,1)+PD.Svecaero(:,1); % BC1 [mmol/m³]
dBC2dt = (Psi*BC2loc' + BC2') + PD.Svec(:,2)+PD.Svecriver(:,2)+PD.Svecaero(:,2); % BC2 [mmol/m³]

dBC1agedt = 1 + (Psi*BC1ageloc' + (PD.Svec(:,1)./BC1loc').*(PD.Svecage(:,1)-BC1ageloc')+... %river+sediment
    (PD.Svecriver(:,1)./BC1loc').*(PD.Svecageriver(:,1)-BC1ageloc')+...
    (PD.Svecaero(:,1)./BC1loc').*(PD.Svecageaero(:,1)-BC1ageloc')); % river DBC Radiocarbon age [d]
dBC2agedt = 1 + (Psi*BC2ageloc' + (PD.Svec(:,2)./BC2loc').*(PD.Svecage(:,2)-BC2ageloc')+... % river+sediment
    (PD.Svecriver(:,2)./BC2loc').*(PD.Svecageriver(:,2)-BC2ageloc')+...
    (PD.Svecaero(:,2)./BC2loc').*(PD.Svecageaero(:,2)-BC2ageloc')); % Radiocarbon age [d] % river+DBC


dwagedt= 1*PD.wagedt'+(Psi*wageloc'); % 0 at surface, +1 everywhere else
dwagedt([1:4])=0;
%% matrix of derivatives

dydtmat = PE.m0;
dydtmat(PE.Ibc1) = dBC1dt; % Ocean Black Carbon 1 [mmolC/m³]
dydtmat(PE.Ibc2) = dBC2dt; % Ocean Black Carbon 2 [mmolC/m³]
dydtmat(PE.Iage1) = dBC1agedt; % Ocean Black Carbon Age 1 [d]
dydtmat(PE.Iage2) = dBC2agedt; % Ocean Black Carbon Age 2 [d]
dydtmat(PE.Iwage) = dwagedt; % water mass age [d]

dydt = dydtmat(PE.Ires)';



