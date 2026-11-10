% Put plug rate data into a table so they can be copied into a csv.
% prepare data for GPR

[plugRate, ] = binofit(numplugs, numtests); 
T7 = table(plugged', atests', repelem(la, length(atests))', VariableNames={'plugged', 'a', 'la'})
Ta = [T1; T2; T3; T4; T5]
Tb = [T6; T7]

%%
writetable(Ta, '~/Documents/GitHub/probabilistic_plugging/PDE rerun plug rate data/PDE_PlugRate_Rerun1.csv');
writetable(Tb, '~/Documents/GitHub/probabilistic_plugging/PDE rerun plug rate data/PDE_PlugRate_Rerun2.csv');


