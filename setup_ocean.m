function Var = setup_ocean(PE, varargin)

% Optional input argument: Circulation strength factor
% Note: for no circulation use 'Psi0' option (different ode file)
Psi_factor  = 0.5; % Default value
if any(strcmp(varargin, 'Psi_factor'))
    ind  = find(strcmp(varargin, 'Psi_factor'));
    Psi_factor = varargin{ind+1};
end

% Default: no one-way flow
oneway = 0;
if any(strcmp(varargin, 'oneway'))
    oneway = 1;
end

%% Ocean parameters

Aoc = PE.Ae*PE.foc; % Ocean area   [m²]
Voc = 1.292e18;     % Ocean volume [m³]
volboxweight = [0.00345,0.006908,0.02072,0.00690,0.1865,0.2714,0.5040];

if PE.nb==7
    % box indices 
    Isfc  = 1:4; % surface box indices
    Idp   = 5:7; % deep box indices
    iAA   = 1; % Antarctic surface
    iSA   = 2; % Subantarctic surface
    iLL   = 3; % Low latitude surface
    iNA   = 4; % N. Atlantic surface
    iTC   = 5; % Thermocline
    iNADW = 6; % North Atlantic Deep Water
    iAABW = 7; % Antarctic Bottom Water
    iacronyms = {'AA','SA','LL','NA','TC','NADW','AABW'};

    % surface boxes
    H(Isfc)  = [250 250 100 250]; % Surface box depths [m]
    A(Isfc)  = [0.05 0.1 0.75 0.1]*Aoc; % Surface box areas [m²]
    V(Isfc)  = A.*H; % Box volume (surface) [m³]
    % thermocline+deep boxes
    H(5) = 900; 
    A(5) = A(3); 
    V(5) = A(5)*H(5);
    A(6:7) = Aoc; %[m²]
    V(6:7) = [0.35 0.65]*(Voc-sum(V)); % [m³]
    H(6:7) = V(6:7)./A(6:7); % edit: itemwise division
%     V = ones(1,7)*(sum(V)/7);
%     warning('volume the same')
    
%     volboxweight = V/sum(V);
    
elseif PE.nb == 2
    Isfc = 1;
    Idp  = 2;
    A = [Aoc Aoc];
    Vsurf = sum([250 250 100 250].*[0.05 0.1 0.75 0.1]*Aoc);
    V = [Vsurf Voc-Vsurf];
%     V = [0.5*Voc 0.5*Voc];
%     warning('same volume for both boxes')
    H = V./A;
    
elseif PE.nb == 1
    Isfc = 1; % surface box
    Idp  = 1;
    A = Aoc;
    V = Voc;
    H = Voc/Aoc;
    
else
    error('Specify box numer nb=1, 2, or 7.')
end

%% circulation (if DoOcnCirc=1)

if PE.nb == 7
    PsiN_o = 20e6; %20*(1e6);  % NADW overturning [Sv = 10^6 m3 s-1]
    PsiS_o = 20e6; %20*(1e6); % AABW overturning [Sv]
    PsiT_o = 100e6; %100*(1e6); % Thermocline overturning [Sv]
    PsiM_o = 40e6;  %40*(1e6); % convective mixing [Sv]

    % PsiN = NADW path [binary, no unit]
    PsiNmat   = zeros(PE.nb,PE.nb);
    PsiNroute = [iNA iNADW iAABW iAA iSA iTC iNA];
    for i = 1:length(PsiNroute)-1; 
        PsiNmat(PsiNroute(i),PsiNroute(i+1)) = 1;
    end

    % PsiS = AABW path
    PsiSmat   = zeros(PE.nb,PE.nb);
    PsiSroute = [iAA iAABW iAA];
    for i = 1:length(PsiSroute)-1; 
        PsiSmat(PsiSroute(i),PsiSroute(i+1))=1;
    end

    % PsiT = path of Thermocline MOC
    PsiTmat   = zeros(PE.nb,PE.nb);
    PsiTroute = [iLL iTC iLL]; 
    for i = 1:length(PsiTroute)-1; 
        PsiTmat(PsiTroute(i),PsiTroute(i+1)) = 1;
    end

    % PsiM = Path of NorthAtlantic Mixing (Two way)
    PsiMmat=zeros(PE.nb,PE.nb);
    PsiMroute=[iNA iNADW iNA];
    for i=1:length(PsiMroute)-1; 
        PsiMmat(PsiMroute(i),PsiMroute(i+1))=1;
    end

    % Total transport matrix (integrating Psi S, N, T and M)
    Psi = PsiTmat'*PsiT_o + PsiSmat'*PsiS_o + PsiNmat'*PsiN_o +...
        PsiMmat'* PsiM_o;
    for i = 1:PE.nb
        Psi(i,i) = -sum(Psi(i,:)); % until here, Psi has units of Sv/10^6 (so not sverdrup, but sverdrup per million!)
    end
    
    Psi_Sv = Psi;
    Psi   = (24*60*60)*Psi./repmat(V',[1 PE.nb]); % divide by volume to get timescale [m³d-1/m³ = d-1]
    Psi_o = Psi*Psi_factor; % by default, factor = 1
    % I think, Psi: column = source box, row = receiving box
    % sum(diag(Psi_o))
    % total circulation in m³/s: sum(-diag(Psi_Sv)) 
    % total default circulation rate in m³/d sum(-diag((24*60*60)*Psi_Sv))
    
elseif PE.nb == 2
    Psi_Downo = 90*(1e6);
    Psi_Upo = 90*(1e6);
    
    % PsiDown = down path [binary, no unit]
    PsiDownmat   = zeros(PE.nb,PE.nb);
    PsiDownroute = [1 2];
    for i = 1:length(PsiDownroute)-1; 
        PsiDownmat(PsiDownroute(i),PsiDownroute(i+1)) = 1;
    end
    
    % PsiUp = Up path [binary, no unit]
    PsiUpmat   = zeros(PE.nb,PE.nb);
    PsiUproute = [2 1];
    for i = 1:length(PsiUproute)-1; 
        PsiUpmat(PsiUproute(i),PsiUproute(i+1)) = 1;
    end

    % Total transport matrix (integrating Psi Down and Up)
    Psi = PsiDownmat'*Psi_Downo + PsiUpmat'*Psi_Upo;
    for i = 1:PE.nb
        Psi(i,i) = -sum(Psi(i,:)); % until here, Psi has units of Sv/10^6 (so not sverdrup, but sverdrup per million!)
    end

    Psi_Sv = Psi;
    Psi   = (24*60*60)*Psi./repmat(V',[1 PE.nb]); % divide by volume to get timescale [m³d-1/m³ = d-1]
    Psi_o = Psi*Psi_factor; % by default, factor = 1
    % I think, Psi: column = source box, row = receiving box
    % sum(diag(Psi_o))
    if oneway == 1 % no backflow from deep into surface box (i.e. oneway flow)
        Psi(:,2) = 0; 
        Psi_o = Psi*Psi_factor;
    end   
    
else
    Psi   = 0; % EDIT: Psi = 1./repmat(V',[1 PE.nb]);
    Psi_o = 0; % EDIT: Psi = 1;
end

%% surface properties

if PE.nb == 7
    T  = [-0.57 5.49 20.42 5.16 10.05 6.52 0.19];%[0 8 22 3]; % surface box temperature (degC)
%     S(Isfc)  = [34 34 36 35.5]; % surface box salinity (psu)
%     TA(Isfc) = [2390 2350 2320 2325]; % Alkalinity [mmol/m3]
elseif PE.nb ==2
    T = [mean([-0.57; 5.49; 20.42; 5.16]) mean([10.05 6.52 0.19])];
%     S(Isfc) = mean([34 34 36 35.5]);
%     TA(Isfc) = mean([2390 2350 2320 2325]);
else
    T = mean([-0.57 5.49 20.42 5.16 10.05 6.52 0.19]);
%     S = mean([34 34 36 35.5]);
%     TA = mean([2390 2350 2320 2325]);
end


%% Compute salinity, temperature and alkalinity
  
% NEWLY COMMENTED
% % preformed matrix
% Apre        = eye(PE.nb,PE.nb);
% Apre(Idp,:) = Psi(Idp,:);

% % PO4
% Amat        = Psi+Qup+Qrem; 
% % as the columns of this matrix would sum to 1, we include another equation
% % indirectly: the mass balance. We know that all boxes' PO4 content must
% % sum to the total initial sum of nutrients. Therefore we make the last row
% % of the matrix contain the fraction of the total volume of the box (to
% % scale them) and make this row sum to pbar, the total global mean of PO4
% % (the other rows of b are 0).
% 
% Amat(end,:) = V/sum(V);
% b      = zeros(1,PE.nb);
% b(end) = Pbar;
% Peq    = (Amat\b')*1e3;
% 
% Ppre, Preg
% bpre       = zeros(PE.nb,1); 
% bpre(Isfc) = Peq(Isfc);
% Ppre       = Apre\bpre;
% Preg       = Peq-Ppre;

% T, S
% bpre       = zeros(PE.nb,1); 
% bpre(Isfc) = T(Isfc);
% T          = (Apre\bpre)';

% NEWLY COMMENTED
% bpre       = zeros(PE.nb,1); 
% bpre(Isfc) = S(Isfc);
% S          = (Apre\bpre)';

% % Age(Isfc) = 0;
% % Age(5)    = V(5)/(PsiT_o+PsiN_o)/PE.spery;
% % Age(6)    = V(6)/(PsiM_o+PsiN_o)/PE.spery;
% % Age(7)    = V(7)/(PsiS_o+PsiN_o)/PE.spery;

%% climate parameters

% cp     = 4185; % specific heat of water [J/kg/K]
% rho_o  = sw_dens(S,T,T*0); % density [kg/m3]
% cm     = cp*rho_o.*(V./A); % specific heat capacity of ocean boxes [J/K/m²?]
% lambda      = 1/(10*PE.spery)*V; 
% lambda(Idp) = 0; % radiative damping 
% dPsidb = 1/3; % scaling between overturning and buoyancy gradient
% 
% % turn off circulation if DoOcnCirc=0
% Psi   = Psi*PS.DoOcnCirc;
% Psi_o = Psi_o*PS.DoOcnCirc;

clear PE

Var=v2struct;



