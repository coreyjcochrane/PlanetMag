% Compare magnetic field measurements from spacecraft against the evaluated
% magnetic field model(s).

% Datasets:
% Voyager 2 MAG: https://doi.org/10.17189/1520034, volume VG2-U-MAG-4-RDR-U1COORDS-48SEC-V1.0

cspice_kclear;
parentName = 'Uranus';
scName = "Voyager 2";
sc = 'Voyager 2';
SEQUENTIAL = 0; % Whether to plot points by index or hours relative to CA

fullOrbFormatSpec = '%24s%5d%2d%2d%1d%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f%3d%8s%[^\n\r]';
disp(['Importing PDS files over ' parentName ' flybys.'])
datFile = fullfile(['MAG/' sc '/vg2_' parentName '_48s_sph.tab']);
orbStr = [parentName ' flyby'];
disp(['Loading ' sc ' MAG data from ' datFile '.'])
fileID = fopen(datFile,'r');
magData = textscan(fileID, fullOrbFormatSpec, inf, 'Delimiter', ',', 'TextType', 'char', 'EndOfLine', '\r\n');
fclose(fileID);

t_UTC = magData{1}';
BrSC = magData{6}';
BthSC = magData{7}';
BphiSC = magData{8}';

LoadSpice(parentName, sc);
Rp_km = 25559;    

nTot = length(t_UTC);
disp(['Converting UTC strings to TDB seconds for all ' num2str(nTot) ' points.'])
ets = cspice_str2et(t_UTC);
disp(['Getting moon and planet distances for all ' num2str(nTot) ' points.'])
[~, rUra_km] = GetMoonDist(sc, parentName, ets);

% Delete measurement times far from Uranus and junk data
BmagSC = sqrt(BrSC.^2 + BthSC.^2 + BphiSC.^2);
moonProx_RP = 0.1;
PlanetMaxDist_RP = 60;
finiteMax_nT = 17e3;
RPunit = ' R_U';
disp(['Excluding all points satisfying at least one of the following:' newline ...
      'Planetocentric distance > ' num2str(PlanetMaxDist_RP) RPunit newline ...
      'Suspect measurements, |B| > ' num2str(finiteMax_nT) 'nT.'])
% Full limits
exclude = find(rUra_km/Rp_km > PlanetMaxDist_RP | BmagSC > finiteMax_nT);
ets(exclude) = [];
BrSC(exclude) = [];
BthSC(exclude) = [];
BphiSC(exclude) = [];

npts = length(ets);
t_h = ets / 3600;
disp(['Getting ' sc ' positions for ' num2str(npts) ' pts.'])
[r_km, theta, phi, xyz_km, spkParent] = GetPosSpice(sc, parentName, t_h);


%% Plot and calculate products
nOpts = 2; nMPopts = 2;
opts = 1:nOpts;
MPopts = 1:(nMPopts + 1); % Add 1 to force noMP model in addition
for opt=opts
    for MPopt=MPopts
        GetBplotAndLsq(ets, t_h, r_km, theta, phi, xyz_km, BrSC, BthSC, BphiSC, ...
            scName, parentName, spkParent, orbStr, opt, MPopt, SEQUENTIAL);
    end
end
