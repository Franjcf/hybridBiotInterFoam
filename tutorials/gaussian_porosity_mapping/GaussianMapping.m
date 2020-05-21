%%Before doing anything: 
	%make sure you have epssBound file (included in this directory)
	%Right now it contains the boundary conditions for the fracturing case,
        % but you can modify it to your liking and the conditions you need 

%% Input parameters and Initializing
clear all

%Input the total number of cells
totalCells = 196256; %49064;
%size of each list
L=totalCells;
    
%initializing the porosity list 
epsslist=zeros(L,1);

%defining the average solid fraction and std
epssAvg = 0.6;
std = 0.05;

%Initalazing Text Files. 
system('cp -r epssBound epss'); %copying dummy files into the files we want

epss = fopen('epss','a+'); %solid fraction file

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
