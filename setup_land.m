
function Var = setup_land(PE,PS)

%% terrestrial parameters

VegName = PS.VegName;

%------------------------------------------------------------------
%---- pools
iLeaf = 1;  % leaf
iRoot = 2;  % root
iWood = 3;  % wood
iLit1 = 4;  % metabolic litter
iLit2 = 5;  % structural litter
iCWDC = 6;  % coarse woody debris
iSOM1 = 7;  % soil organic matter 1
iSOM2 = 8;  % soil organic matter 2
iSOM3 = 9;  % soil organic matter 3;


%------------------------------------------------------------------

Ala = PE.Ae*PE.fla; % Land area [m²]


%------------------------------------------------------------------
switch VegName  % This loop will switch between biomes and global mean
    
    
%---------------------Global Mean Case ---------------------------%
  case 'Global'
    Q10_resp  = 2;  
    beta_550  = 1.25;
    beta_fert = ( beta_550 - 1 ) / log( 550 / 365 ); 
    NPP_o = 50;    % Net Primary Production [Pg C/yr]
    NPP_o = NPP_o*1e15/12/PE.spery/Ala; % convert to molC/m2/s
    Rh_o = NPP_o; % steady state (preindustrial)



    %------------------------------------------------------------------
    %---- NPP partitioning to pools (leaf, root, wood) [fraction]
    bcoef    = zeros(PE.np,1);
    bcoef(iLeaf) = 0.30;      % leaf
    bcoef(iRoot) = 0.50;      % root
    bcoef(iWood) = 0.20;      % wood

    %------------------------------------------------------------------
    %---- base turnover rate for pool j (tau is in years)

    tauL(iLeaf)=1;
    tauL(iRoot)=10;
    tauL(iWood)=40;
    tauL(iLit1)=0.5;
    tauL(iLit2)=0.5;
    tauL(iCWDC)=50;
    tauL(iSOM1)=0.5;
    tauL(iSOM2)=2.5;
    tauL(iSOM3)=303;

%---------------------Tropical Forest Case ---------------------------%
  case 'TRF'
    Q10_resp  = 2;  
    beta_550  = 1.25;
    beta_fert = ( beta_550 - 1 ) / log( 550 / 365 ); 
    NPP_o = 200;    % Net Primary Production [g C/m2/yr]
    NPP_o = NPP_o/12/PE.spery; % convert to molC/m2/s
    Rh_o = NPP_o; % steady state (preindustrial)

    %------------------------------------------------------------------
    %---- NPP partitioning to pools (leaf, root, wood) [fraction]
    bcoef    = zeros(PE.np,1);
    bcoef(iLeaf) = 0.30;      % leaf
    bcoef(iRoot) = 0.30;      % root
    bcoef(iWood) = 0.40;      % wood

    %------------------------------------------------------------------
    %---- base turnover rate for pool j (tau is in years)

    tauL(iLeaf)=2;
    tauL(iRoot)=2;
    tauL(iWood)=100;
    tauL(iLit1)=0.5;
    tauL(iLit2)=0.5;
    tauL(iCWDC)=10;
    tauL(iSOM1)=0.5;
    tauL(iSOM2)=2.5;
    tauL(iSOM3)=50;
        
%---------------------Temperate Forest Case --------------------------%    
    case 'TempForest'
               
    Q10_resp  = 2;  
    beta_550  = 1.25;
    beta_fert = ( beta_550 - 1 ) / log( 550 / 365 ); 
    NPP_o = 80;    % Net Primary Production [g C/m2/yr]
    NPP_o = NPP_o/12/PE.spery; % convert to molC/m2/s
    Rh_o = NPP_o; % steady state (preindustrial)

    %------------------------------------------------------------------
    %---- NPP partitioning to pools (leaf, root, wood) [fraction]
    bcoef    = zeros(PE.np,1);
    bcoef(iLeaf) = 0.30;      % leaf
    bcoef(iRoot) = 0.30;      % root
    bcoef(iWood) = 0.40;      % wood

    %------------------------------------------------------------------
    %---- base turnover rate for pool j (tau is in years)

    tauL(iLeaf)=2;
    tauL(iRoot)=2;
    tauL(iWood)=30;
    tauL(iLit1)=0.5;
    tauL(iLit2)=0.5;
    tauL(iCWDC)=25;
    tauL(iSOM1)=0.5;
    tauL(iSOM2)=2.5;
    tauL(iSOM3)=300;
        
%---------------------Boreal Forest Case -----------------------------%
    case 'BorealForest'
        
    Q10_resp  = 2;  
    beta_550  = 1.25;
    beta_fert = ( beta_550 - 1 ) / log( 550 / 365 ); 
    NPP_o = 35;    % Net Primary Production [g C/m2/yr]
    NPP_o = NPP_o/12/PE.spery; % convert to molC/m2/s
    Rh_o = NPP_o; % steady state (preindustrial)

    %------------------------------------------------------------------
    %---- NPP partitioning to pools (leaf, root, wood) [fraction]
    bcoef    = zeros(PE.np,1);
    bcoef(iLeaf) = 0.30;      % leaf
    bcoef(iRoot) = 0.30;      % root
    bcoef(iWood) = 0.40;      % wood

    %------------------------------------------------------------------
    %---- base turnover rate for pool j (tau is in years)

    tauL(iLeaf)=2;
    tauL(iRoot)=2;
    tauL(iWood)=35;
    tauL(iLit1)=0.5;
    tauL(iLit2)=5;
    tauL(iCWDC)=50;
    tauL(iSOM1)=0.5;
    tauL(iSOM2)=10;
    tauL(iSOM3)=500;
        
%---------------------Grassland Case ---------------------------------%
    case 'Grass'
        
    Q10_resp  = 2;  
    beta_550  = 1.25;
    beta_fert = ( beta_550 - 1 ) / log( 550 / 365 ); 
    NPP_o = 20;    % Net Primary Production [g C/m2/yr]
    NPP_o = NPP_o/12/PE.spery; % convert to molC/m2/s
    Rh_o = NPP_o; % steady state (preindustrial)

    %------------------------------------------------------------------
    %---- NPP partitioning to pools (leaf, root, wood) [fraction]
    bcoef    = zeros(PE.np,1);
    bcoef(iLeaf) = 0.50;      % leaf
    bcoef(iRoot) = 0.50;      % root
    bcoef(iWood) = 0.0;      % wood

    %------------------------------------------------------------------
    %---- base turnover rate for pool j (tau is in years)

    tauL(iLeaf)=2;
    tauL(iRoot)=2;
    tauL(iWood)=1;
    tauL(iLit1)=0.5;
    tauL(iLit2)=0.5;
    tauL(iCWDC)=1;
    tauL(iSOM1)=0.5;
    tauL(iSOM2)=2.5;
    tauL(iSOM3)=10;
        
end  



%------------------------------------------------------------------
%- the remaning variables will be the same for all biomes
%------------------------------------------------------------------
kbase      = zeros(PE.np,1);
kbase(iLeaf) = 1 / tauL(iLeaf);
kbase(iRoot) = 1 / tauL(iRoot);
kbase(iWood) = 1 / tauL(iWood);

kbase(iLit1) = 1 / tauL(iLit1);
kbase(iLit2) = 1 / tauL(iLit2);
kbase(iCWDC) = 1 / tauL(iCWDC);

kbase(iSOM1) = 1 / tauL(iSOM1);
kbase(iSOM2) = 1 / tauL(iSOM2);
kbase(iSOM3) = 1 / tauL(iSOM3);

kbase = kbase/PE.spery; % convert to 1/sec
krate = diag(kbase);    % scaled turnover rate

%------------------------------------------------------------------
%---- fractional carbon flow from pool j to pool i (pathf(i,j))
pathf      = zeros(PE.np,PE.np);
pathf(iLit1,iLeaf) = 0.6;
pathf(iLit2,iLeaf) = 0.4;

pathf(iLit1,iRoot) = 0.6;
pathf(iLit2,iRoot) = 0.4;

pathf(iCWDC,iWood) = 1.0;

pathf(iSOM1,iLit1) = 1.0;

pathf(iSOM2,iLit2) = 0.15;
pathf(iSOM1,iLit2) = 1. - pathf(iSOM2,iLit2);

pathf(iSOM2,iCWDC) = 0.25;
pathf(iSOM1,iCWDC) = 1. - pathf(iSOM2,iCWDC);

pathf(iSOM3,iSOM1) = 0.01;
pathf(iSOM2,iSOM1) = 1. - pathf(iSOM3,iSOM1);

pathf(iSOM3,iSOM2) = 0.005;
pathf(iSOM1,iSOM2) = 1. - pathf(iSOM3,iSOM2);

pathf(iSOM1,iSOM3) = 1.0;

%------------------------------------------------------------------
%---- fractional respiration loss for carbon flow from pool j to pool i
respf      = zeros(PE.np,PE.np);

respf(iLit1,iLeaf) = 0.;
respf(iLit2,iLeaf) = 0.;

respf(iLit1,iRoot) = 0.;
respf(iLit2,iRoot) = 0.;

respf(iCWDC,iWood) = 0.;

respf(iSOM1,iLit1) = 0.55;

respf(iSOM2,iLit2) = 0.30;
respf(iSOM1,iLit2) = 0.55;

respf(iSOM2,iCWDC) = 0.30;
respf(iSOM1,iCWDC) = 0.45;

respf(iSOM3,iSOM1) = 0.;
respf(iSOM2,iSOM1) = 0.45;

respf(iSOM3,iSOM2) = 0.;
respf(iSOM1,iSOM2) = 0.55;

respf(iSOM1,iSOM3) = 0.55;

%------------------------------------------------------------------
%---- fractional carbon flow from pool j that enters pool i
acoef = diag(-1.0 .* ones(PE.np,1)) + pathf .* (1. - respf);

%% save all variable into structure

clear PE

Var=v2struct;




