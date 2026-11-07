% Put plug rate data into a table so they can be copied into a csv.
% prepare data for GPR

[plugRate, ] = binofit(numplugs, numtests); 
T5 = table(plugRate', a_values', repelem(la, length(plugRate))', VariableNames={'plugRate', 'a', 'la'})
T0 = [T1; T2; T3; T4; T5]