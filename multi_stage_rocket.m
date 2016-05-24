%% ***** Orbital Maneuver Calculations *******
clear all

% Choose Destination
%********************

disp('Choose destination, type the corresponding number.')
disp('1. Mun')
disp('2. Duna')
promt = 'Destination = ';
des = input(promt);

% Constants
g0 = 9.82; %(m/s^2) Gravitational acceleration constant
mu = 3.5316e3; % (Kerbin) | 3.986012e5 (Earth); % Gravitational Parameter (km^3/s^2)
r1 = 700; % 6697.0523 (Earth); % Parking Orbit Radius, LEO (km)

%% ***** Trans Mun Insertion *********

if des == 1
disp('Choose mission type, type the corresponding number.')
disp('1. Munar Orbit')
disp('2. Circummunar Flight')
prompt = 'Mission = ';

mish = input(prompt);
    
% Initial Conditions
phi1 = 0; % Flight Path Angle (deg)
lampp = 30; % Patch Angle (deg)
v1 = 3.09268; % 10.8462 (Earth); % Transfer Orbit departure Velocity (km/s)

% Constants
Rs = 2430; % (Mun) | 66300 (Moon); % Sphere of Influence Radius (km)
D = 12000; % (Mun) | 384400 (Moon); % Mean Munar Orbit Radius (km)
rm = 200; % (Mun) | 1738 (Moon); % Radius of Mun (km)
vm = 0.5425; % (Mun) | 1.018 (Moon); % Velocity of Mun relative to center of Kerbin (km/s)
wm = vm/rm; % (Mun) | 1.5178e-4 (Moon); % Mun Angular Velocity (deg/s)
mum = 65.138398; % (Mun) | mu/81.3 (Moon); % Mun Gravitational Parameter (km^3/s^2)
end

% Geocentric Departure Orbit
%****************************

% Orbit Departure
vcs1 = sqrt(mu/r1) % Lower Earth Orbit Velocity (km/s)
disp('km/s')
dv1 = v1 - vcs1 % Departure Delta-V (km/s)
disp('km/s')
dv1_m = dv1*1000; % Departure Delta-V (m/s)
Emt = (v1^2)/2 - mu/r1; % Transfer Orbit Energy (km^2/s^2)
%Emt = -0.011
Hm1 = r1*v1; % Angular Momentum (km^3/s)

% Orbit Patch Point
rpp = sqrt(Rs^2 + D^2 - 2*Rs*D*cosd(lampp)) % Patch Point Radius (km)
disp('km')
vpp = sqrt(2*(Emt + mu/rpp)) % Patch Point Velocity (km/s)
disp('km/s')
phipp = acosd(Hm1/(rpp*vpp)) % Patch Point Phase Angle (deg)
disp('deg')
gampp = asind((Rs/rpp)*sind(lampp)) % SOI Arrival Angle (deg)
disp('deg')

% Conditions at Patch Point
r2 = Rs; % Selenocentric Radius (km)
v2 = sqrt(vpp^2 + vm^2 - 2*vpp*vm*cosd(phipp-gampp)); % Selenocentric Velocity (km/s)
% v2_km_s = v2*(6378.145/806.8118744)
ep2 = asind((vm/v2)*cosd(lampp)-(vpp/v2)*cosd(lampp+gampp-phipp)); % Selenocentric Velocity Angle (deg)

% Selenocentric Arrival Orbit
%*****************************

% Orbital Elements
Emm = (v2^2)/2 - mum/r2; % Munar Orbit Energy (km^2/s^2)
Hmm = r2*v2*sind(ep2); %  Munar Orbit Angular Momentum (km^3/s)
pm = (Hmm^2)/mum; % (km)
em = sqrt(1 + 2*((Emm*(Hmm^2))/(mum^2))); % Munar Orbit Eccentricity

% Conditions at Periselenium
rp = pm/(1+em) % Periselenium Radius (km)
disp('km')
% rp_km = rp*6378.145
vp = sqrt(2*(Emm + mum/rp)) % Periselenium Velocity (km/s)
disp('km/s')
hp = rp - rm % Minimum Distance of Approach (km)
disp('km')
% hp_km = rp_km - rm*6378.145

% Mission Objective
if mish == 1
    vcs_mun = sqrt(mum/rp) % Circular Munar Orbit Velocity at Periselenium (km/s)
    disp('km/s')
    dv2 = vp - vcs_mun % Circular Munar Orbit Insertion Delta-V (km/s)
    disp('km/s')
    dv2_m = dv2*1000; % Insertion Delta-V (m/s)
elseif mish == 2
    
end

%% Calculate Time of Flight and Phase Angle

% Orbital Eliments
pt = (Hm1^2)/mu; % (km)
at = -mu/(2*Emt); % Semi-major Axis (km)
et = sqrt(1 - pt/at); % Eccentricity 

% Calculate Time of Flight
vA1 = 0; % True Anomaly, Departure (deg)
vApp = acosd((pt - rpp)/(rpp*et)); % True Anomaly, Patch Point (deg)
% disp('deg')
EA1 = acosd((et + cosd(vA1))/(1 + et*cosd(vA1))); % Eccentric Anomaly, Dep (deg)
EA1_rad = EA1*(pi/180); % (rad)
EApp = acosd((et + cosd(vApp))/(1 + et*cosd(vApp))); % Eccentric Anom, PP (deg)
EApp_rad = EApp*(pi/180); % (rad)
tof = sqrt((at^3)/mu)*((EApp_rad-et*sin(EApp_rad))-(EA1_rad-et*sin(EA1_rad))) % Time of Flight (s)
disp('sec')
tof_hrs = tof*(1/3600) % (hrs)
disp('hrs')
gam1 = vApp - vA1 - gampp - wm*tof % Phase Angle, Departure (deg)
disp('deg')


%% ***** Rocket Design ********************

% n = 2; % No. of stages

% Isp1 = 280; %(sec) Specific Impulse
% Isp2 = 300;

Isp1 = 280; %(sec) Specific Impulse
Isp2 = 320;

Vex = [Isp1*g0,Isp2*g0]; % Exhaust velocity (m/s)

dVe_LEO = 4700; %(m/s) Delta-V from Launch to 100 km altitude orbit

if des == 1
    if mish == 1
        dVe = dVe_LEO + dv1_m + dv2_m % Total Delta-V to mission objective (m/s)
        disp('m/s')
        n = 2; % No. of stages
    elseif mish == 2
    
    end
else
    n = 2; % No. of stages
end

%*** Example ***
% Vex = [3528,3528,3528];
% dVe = 9077;
% n = 3;
%***************

% % Vorb = 2245.5; %(m/s) Required orbital velocity
% % VrE = 465; %(m/s) Speed of rotation of Earth at Equator (427 at Kennedy Space Center)
% % Vgl = 1200; %(m/s) Gravitational velocity losses
% % Va = 1500; %(m/s) Aerodynamic velocity losses
% % 
% % Vn = Vorb - VrE + Vgl + Va; %(m/s) Final velocity

% Masses in Tons
% MLt = 1.5; %(tons) Payload mass
MLt = 6.1;
% M0t = [30.7,4.5]; %(tons) Total initial mass of the ith stage
% Mst = [6.2,1.6]; %(tons) Structural mass of the ith stage
% Mpt = [20.0,2.0]; %(tons) Propellant mass of the ith stage

tkg = 1000; % Tonnes to Kg conversion factor

% Convert to kg
ML = MLt*tkg; %(kg) Payload mass
% M0 = [M0t(1)*tkg,M0t(2)*tkg,ML] %(kg) Total initial mass of the ith stage
% Ms = [Mst(1)*tkg,Mst(2)*tkg] %(kg) Structural mass of the ith stage
% Mp = [Mpt(1)*tkg,Mpt(2)*tkg] %(kg) Mass of propellant in the ith stage

%% ***** Structural Coefficients *********

% Fuel Tank Masses
%******************

% Small
%-------

% Short                     Long
MoS_s = 2.4*tkg;            MoL_s = 5.8*tkg;            
MpS_s = 1.0*tkg;            MpL_s = 4.0*tkg;            
MsS_s = MoS_s - MpS_s;      MsL_s = MoL_s - MpL_s;      

s_eff_tank(1,1) = MsS_s/MoS_s;
s_eff_tank(1,2) = MsL_s/MoL_s;

% Medium
%--------

% Short                     Long
MoS_m = 15.0*tkg;           MoL_m = 42.0*tkg;   
MpS_m = 8.0*tkg;            MpL_m = 32.0*tkg;    
MsS_m = MoS_m - MpS_m;      MsL_m = MoL_m - MpL_m;      

s_eff_tank(2,1) = MsS_m/MoS_m;
s_eff_tank(2,2) = MsL_m/MoL_m;

% Large
%-------

% Short                     Long
MoS_l = 35.3*tkg;           MoL_l = 96.0*tkg;           
MpS_l = 18.0*tkg;           MpL_l = 72.0*tkg;           
MsS_l = MoS_l - MpS_l;      MsL_l = MoL_l - MpL_l;      

s_eff_tank(3,1) = MsS_l/MoS_l;
s_eff_tank(3,2) = MsL_l/MoL_l;

%*** Example ***
% s_eff = [0.1,0.1,0.1];
%***************
% s_eff = [0.05,0.07]; % Structural coefficient
% s_eff = [0.25,0.43];
% s_eff = [0.2381,0.4667]; % Small Size
  s_eff = [0.2130,0.2689]; % Medium Size
%         Stage 1  Stage 2
% s_eff = [0.2130,0.2638]; % Found fromt test-1-mun
 

% s_eff = [s_eff_tank(1,2),s_eff_tank(1,1)] % Data order is flipped




%% ***** Rocket Assembly Calculations **********

for i = 1:n
%    s_eff(i) = Ms(i)/(M0(i) - M0(i+1))
%    s_eff(i) = Ms(i)/(Mp(i) + Ms(i))
    R(i) = exp(dVe/(n*Vex(i))); % Mass ratio
    a(i) = Vex(i)*(1 - s_eff(i)*R(i)); % Lagrange multiplier (m/s)
    lam(i) = (1 - s_eff(i)*R(i))/(R(i) - 1); % Payload ratio
    gam(i) = ((1 - s_eff(i)*R(i))/((1 - s_eff(i))*R(i)))^n; % Payload fraction for each stage
end

% A = a
% Lam = lam
% Mass_R = R
% Gam_i = gam
Gam = sum(gam) % Overall payload fraction


% Determine Total Masses
M0(1) = ML/Gam;
M0(2) = M0(1)*(lam(1)/(1 + lam(1)));
%M0(3) = M0(2)*(lam(2)/(1 + lam(2)));
M0(3) = MLt*tkg;

% Determine Propellent & Structural Masses
for i = 1:n
    % Mp2(i) = M0(i) - M0(i)/R(i);
    % Ms(i) = (s_eff(i)*Mp(i))/(1 - s_eff(i));
    Ms(i) = s_eff(i)*(M0(i) - M0(i+1));
    Mp(i) = Ms(i)/s_eff(i) - Ms(i);
end

MLt
out_Mst = Ms/tkg
disp('tons')
out_Mpt = Mp/tkg
disp('tons')
out_M0t = M0/tkg
disp('tons')
