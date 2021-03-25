% read the data

%a = [180 154 143 134 127 120 114 107 102 96 90 84 78 73 66 60 53 46 37 26 0];
En = [60 55 50 45 40 35 30 25 20 15 10 5];
a = [180 154 143 134 127 120 114 107 102 96];


for k =1:12;
for n =1:10;
    %close all; 
    clear particle; 
    filename = ['/media/andrea/FREECOM HDD/Orbits/ORBIT_29880_at_0.255_s_' num2str(En(k)) '_keV_' num2str(a(n)) '_deg_60_deg_R0.90_Z0.00.mat'];
   %filename = ['/home/andrea/Documents/MASTORBIT/EPS2017/pre_2ndsawtooth/pre_2ndsawtooth_0.9/ORBIT_29880_at_0.255_s_' num2str(En(k)) '_keV_' num2str(a(n)) '_deg_R0.90_Z0.00.mat'] 
    load(filename);
    fprintf('%d keV and %d deg\n', En(k), a(n));
    MOrbit_plot_selection(particle, EQLST, 1);
    
    fprintf('done\n'); fflush(stdout());
    input('Next ?', 'answer');
endfor
endfor