
function Var = setup_run(PE, varargin)

% Optional input argument: Number of boxes
nb  = 7; % Default value
if any(strcmp(varargin, 'nb'))
    ind  = find(strcmp(varargin, 'nb'));
    nb = varargin{ind+1};
end

% Optional input argument: Length of spin up period [y]
yspin  = 30000; % Default value 400
if any(strcmp(varargin, 'yspin'))
    ind  = find(strcmp(varargin, 'yspin'));
    yspin = varargin{ind+1};
end

%% Domains and processes to turn on or off (1="on", 0="off")

% mean state switches - to run land only, set ocean values to 0

runName = 'uncoupled';
% 
% DoOcn = 1; % include ocean 
% DoOcnBio  = 0%1; % include ocean biological pump (uptake/remin of Corg)
% DoOcnSol  = 0%1; % include ocean solubility pump (SST variation)
% DoOcnCirc = 0%1; % include ocean solubility pump (Circulation variation)
% 
% % carbon-climate coupling (CCC) processes
% DoRadCO2 = 0%1; % let climate respond to CO2 changes
% CCC_OT   = 0%1; % allow CCC via temperature on ocean C (e.g. solubility, circulation)
% CCC_OC   = 0%1; % allow CCC via atm CO2 effect on ocean C (e.g. disequilibrium, buffer factor)

%% emissions scheme

% specify co2 scenario
escheme = 'zero'; % no emission
% escheme = 'rcp85'; % emissions trajectory % EDIT
%escheme = 'pulse'; % emissions trajectory [WARNING: NOT YET CODED]
%escheme = 'ramp1'; % co2 trajectory [WARNING: NOT YET CODED]
cmax = 4; % co2 change factor used in ramp scenarios (e.g. 4 means 4x[patm0])

%% climate forcing 

% specify climate forcing 
RFamp = 0; % amplitude of imposed RF oscillations [W/m2]
RFper = 100; % period of oscillatory forcing [yr]

%% Terrestrial ecosystem and forcing

%---->>>>> Specify Ecosystem type here! <<<<----%
% Options: 'Global', 'TRF', 'TempForest','BorealForest','Grass'
VegName='Global';

% specify NPP forcing 
NPPamp = 0; % fractional amplitude of imposed NPP oscillations (0.5 means 50% change)
NPPper = 4; % period of oscillatory forcing [yr]


%% time parameters

% the ocean is first run into equilibrium with the atmosphere prescribed a
% certain carbon content, during this time the atmosphere will keep this
% prescribed value, no matter what the ocean does. So during this time, it
% is not carbon-balanced.
% ypert is the time when the emission scheme starts and the atmosphere gets
% coupled to the ocean. so from this time point on, there should be a mass
% balance. 

% yspin = 5000; % length of spinup period in years (needed for ocean)
yrun  = 100; % time at end (calendar year) [no influence for escheme zero]
yend  = yspin + yrun;
% assert(yend>ypert, 'Calender year at end must be > year of pertubation')
ytot=yspin+yrun; % changed by Sinikka 07/2021
%% radiative and NPP forcing 

% forcing is external (non-CO2, non-Temp mechanisms)
% may be stochastic or periodic (requires manual switch)

Yint = yspin : 0.1 : yend; % year of perturbation time series

%RFint = RFamp * (rand(length(Yint),1)-0.5)*2;
RFint = RFamp * (sin(((Yint-yspin)*2*pi/(RFper))));

NPPint = NPPamp * (sin(((Yint-yspin)*2*pi/(NPPper))));

%% save all variable into structure

clear PE

Var=v2struct;


