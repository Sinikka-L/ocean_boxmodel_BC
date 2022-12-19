function [pco2local,pHlocal,fflocal] = calc_pco2(t,s,ta,c,phg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the pCO2 of seawater
% assuming constant PO4 (0uM) and SiO3 (40uM)
% INPUTS:
%   t -> temperature [deg C]
%   s -> salinity [psu]
%   ta -> total alkalinity [uM]
%   c -> DIC [uM]
%   phg -> guess at pH solution [unitless]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pt  = 0e-3;
sit = 40.0e-3;
tk  = 273.15 + t;
tk100 = tk/100.0;
tk1002= tk100.*tk100;
invtk = 1.0./tk;
dlogtk= log(tk);
is    = 19.924*s./(1000.-1.005*s);
is2   =is.*is;
sqrtis=sqrt(is);
s2    =s.*s;
sqrts =sqrt(s);
s15   =s.^1.5;
scl   =s/1.80655;

fflocal = exp(-162.8301 + 218.2968./tk100  + ...
     90.9241*log(tk100) - 1.47696*tk1002 + ...
     s .* (.025695 - .025225*tk100 + ...
	  0.0049867*tk1002));

k0local = exp(93.4517./tk100 - 60.2409 + ...
     23.3585 * log(tk100) + ...
     s .* (0.023517 - 0.023656*tk100 + ...
     0.0047036*tk1002));
 
k1local = 10.^(-1*(3670.7*invtk - ...
     62.008 + 9.7944*dlogtk - ...
     0.0118 * s + 0.000116*s2));

k2local = 10.^(-1*(1394.7*invtk + 4.777 - ...
     0.0184*s + 0.000118*s2));

kblocal = exp((-8966.90 - 2890.53*sqrts - 77.942*s + ...
     1.728*s15 - 0.0996*s2).*invtk + ...
     (148.0248 + 137.1942*sqrts + 1.62142*s) + ...
     (-24.4344 - 25.085*sqrts - 0.2474*s) .* ...
     dlogtk + 0.053105*sqrts.*tk);

k1plocal = exp(-4576.752*invtk + 115.525 - ... 
     18.453*dlogtk + ...
     (-106.736*invtk + 0.69171).*sqrts + ...
	       (-0.65643*invtk - 0.01844).*s);

k2plocal = exp(-8814.715*invtk + 172.0883 - ...
     27.927*dlogtk + ... 
     (-160.340*invtk + 1.3566) .* sqrts + ... 
     (0.37335*invtk - 0.05778) .* s);

k3plocal = exp(-3070.75*invtk - 18.141 + ... 
     (17.27039*invtk + 2.81197) .* ... 
     sqrts + (-44.99486*invtk - 0.09984) .* s);

ksilocal = exp(-8904.2*invtk + 117.385 - ...
     19.334*dlogtk + ...
     (-458.79*invtk + 3.5913) .* sqrtis + ...
     (188.74*invtk - 1.5998) .* is + ... 
     (-12.1652*invtk + 0.07871) .* is2 + ... 
	       log(1.0-0.001005*s));

kwlocal = exp(-13847.26*invtk + 148.9652 - ...
     23.6521*dlogtk + ...
     (118.67*invtk - 5.977 + 1.0495 * dlogtk) .* ... 
	      sqrts - 0.01615 * s);

kslocal = exp(-4276.1*invtk + 141.328 - ...
     23.093*dlogtk + ...
     (-13856*invtk + 324.57 - 47.986*dlogtk).*sqrtis + ...
     (35474*invtk - 771.54 + 114.723*dlogtk).*is - ...
     2698*invtk.*is.^1.5 + 1776*invtk.*is2 + ...
	      log(1.0 - 0.001005*s));

kflocal = exp(1590.2*invtk - 12.641 + 1.525*sqrtis + ...
     log(1.0 - 0.001005*s) + ...
	      log(1.0 + (0.1400/96.062)*(scl)./kslocal) );

btlocal = 0.000232 * scl/10.811;
stlocal = 0.14 * scl/96.062;
ftlocal = 0.000067 * scl/18.9984;

pHlocal = phg;
permil=1.0/1024.5;
pt=pt*permil;
sit=sit*permil;
ta=ta*permil;
c = c*permil;

%%%%%%%%%%%%%%%%%%%%%
%% start iteration %%
%%%%%%%%%%%%%%%%%%%%%

phguess = pHlocal;
hguess = 10.0.^(-phguess);
bohg = btlocal.*kblocal./(hguess+kblocal);
stuff = hguess.*hguess.*hguess ...
     + (k1plocal.*hguess.*hguess) ...
     + (k1plocal.*k2plocal.*hguess) ...
     + (k1plocal.*k2plocal.*k3plocal);
h3po4g = (pt.*hguess.*hguess.*hguess) ./ stuff;
h2po4g = (pt.*k1plocal.*hguess.*hguess) ./ stuff;
hpo4g  = (pt.*k1plocal.*k2plocal.*hguess) ./ stuff;
po4g   = (pt.*k1plocal.*k2plocal.*k3plocal) ./ stuff;

siooh3g = sit.*ksilocal ./ (ksilocal + hguess);

cag = ta - bohg - (kwlocal./hguess) + hguess ... 
     - hpo4g - 2.0*po4g + h3po4g ...
     - siooh3g;

gamm  = c./cag;
hnew = 0.5*(-k1local.*(1-gamm)+sqrt((k1local.^2).*(1 - gamm).^2 ... 
   +4*k1local.*k2local .*(2*gamm - 1) ) );

pHlocal_new = -log10(hnew);
pHlocal = pHlocal_new;

pco2local = c./fflocal./(1.0 + (k1local./hnew) + ... 
   (k1local.*k2local./(hnew.*hnew)) );
fflocal = fflocal/permil;
