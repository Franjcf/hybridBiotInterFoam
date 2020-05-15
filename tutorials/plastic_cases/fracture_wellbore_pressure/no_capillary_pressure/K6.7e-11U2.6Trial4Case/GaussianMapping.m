%%Before doing anything:

%place all images in matlab directory
%place a k and a epsf into directory
%do pre-image analysis and look for the appropiate thresholds to set 

%% Input parameters and Initializing
clear all

%Input the total number of cells
totalCells = 196256; %49064;
%size of each list
L=totalCells;
    
%initializing the porosity list 
epsslist=zeros(L,1);

%defining the average and std
epssAvg = 0.4;%    (fracture verification)
std = 0.05;% 0.075 (fracture verification)

%epssAvg = 0.1036; (fracture verification)
%std = 0.01; (fracture verification

%Initalazing Text Files. 
system('cp -r epssBound epss'); %copying dummy files into the files we want

epss = fopen('epss','a+'); %porosity file

fprintf(epss,'%15.0f\n',L); %total number of values at top

fprintf(epss,'%s\n','(');

%%making random numbers
rng('shuffle');

%% Making the Lists    
for i = 1:L %number of cells
    epsslist(i)= min(0.99,abs(normrnd(epssAvg,std)));
end

%% Writing text files
fprintf(epss,'%1.5e\n',epsslist);
fprintf(epss,'%s\n',');');
fprintf(epss,'%s\n','// ************************************************************************* //');

%closing file
fclose(epss);
fclose('all');

%% print usefull stuff
fprintf('Total number of data points: %.0f\n',L);
