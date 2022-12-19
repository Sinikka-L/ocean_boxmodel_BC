function Var = setup_model(varargin)

nd = 5; % number state variable 'domains' (Oc DOM, Oc Bacteria, Oc DOM age)
% % nd = 4; % number state variable 'domains' (Oc DOM, Oc Bacteria, Oc DOM age)

% Optional input argument: Number of boxes
nb  = 7; % Default value
if any(strcmp(varargin, 'nb'))
    ind  = find(strcmp(varargin, 'nb'));
    nb = varargin{ind+1};
end

%added by Sinikka to make NPP in setup_land work
np=1;
%% Earth parameters 

re       = 6371e3; % Earth radius [m]
Ae       = 4*pi*re^2; % Earth area [m²]
fla      = 0.3; % fraction land area
foc      = 1-fla; % fraction ocean area
ps       = 1013.5 * (1e2); % mean surface pressure [Pa]

%% Structure and Indices of state variable arrays

% generic matrices
m0 = zeros(nb,nd)*nan; 
m1 = m0+1;

% indices of pools in generic matrices
%Idom = sub2ind(size(m0),1:nb,zeros(1,nb)+1);
%Ibac = sub2ind(size(m0),1:nb,zeros(1,nb)+2);
%Iage = sub2ind(size(m0),1:nb,zeros(1,nb)+3); 
%Ires = cat(2,Idom,Ibac,Iage);
% % Ierr = sub2ind(size(m0),1,4);
% % Ierr = sub2ind(size(m0),1:nb,zeros(1,nb)+4); 
% % Ires = cat(2,Idom,Ibac,Iage,Ierr);
Ibc1= sub2ind(size(m0),1:nb,zeros(1,nb)+1);
Ibc2= sub2ind(size(m0),1:nb,zeros(1,nb)+2);
Iage1= sub2ind(size(m0),1:nb,zeros(1,nb)+3); 
Iage2= sub2ind(size(m0),1:nb,zeros(1,nb)+4); 
Iwage= sub2ind(size(m0),1:nb,zeros(1,nb)+5); 
Ires = cat(2,Ibc1,Ibc2,Iage1,Iage2,Iwage);

% indices of pools in vector of m0(Ires)
Jbc1 = 0*nb+1:1*nb; 
Jbc2 = 1*nb+1:2*nb;
Jage1 = 2*nb+1:3*nb;
Jage2 = 3*nb+1:4*nb;
Jwage = 4*nb+1:5*nb;
% % Jerr = 3*nb+1:4*nb;

spery=24*3600*365; % econds per year, added by Sinikka 07/2021
%% save all variable into structure

Var=v2struct;




