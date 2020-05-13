%%To run this on bash type the following, makes nure you have epssBound in te same directory.
% This will create an solid fraction profile that you can then placein the 0/ directory.
%makse sure you have matlab module loaded. Otherwise just open script on matlab and run. 

%type this: matlab -nodisplay -nosplash -nodesktop -r "run('GaussianMapping.m');exit;" 

%% Input parameters and Initializing
clear all

%Input the total number of cells
totalCells = 150000; %49064;
%size of each list
L=totalCells;
    
%initializing the porosity list 
epsslist=zeros(L,1);

%defining the average and std
epssAvg = 0.6;%    (fracture verification)
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
