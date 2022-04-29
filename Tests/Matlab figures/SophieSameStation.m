
%Time between observation in seconds
t=[60,50,40];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%SOPHIE NUSAT-7 SAME LOCATION%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Latitude at observation time
phi=[-54.371917,-12.04181,1.90617,-53.92713,14.64296,20.35315,37.7887];

%earth radius at latitude: https://rechneronline.de/earth-radius/
R=[6364.049,6377.214,6378.114,6364.207,6376.781,6375.569,6370.15]


%Measured height, i.e. r2-earth radius at lattitude
h_computed=[];

%Measured velocity
v_computed=[];

%Measured inclination
i_computed=[97.1997];

%Measured right ascension
Omega_computed=[182.2854];

%Measured eccentricity
e_computed=[1.0062014];

%Measured mean motion
n_computed=[15.17766642];

%TLE: https://orbit.ing-now.com/satellite/45017/2020-003b/nusat-7/
%Comparison TLE 5h old
%Reference value
h_ref=[482.33,478.95,476.14,492.13];
v_ref=[7.6248,7.6205,7.6231,7.6138];
i_ref=97.2200;
Omega_ref=182.1223;
e_ref=0.0012960;
n_ref=15.33717004127678;

