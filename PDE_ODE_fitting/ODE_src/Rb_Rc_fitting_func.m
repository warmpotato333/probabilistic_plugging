% Reed Ogrosky, Aya Yu
% Created March 2, 2026
% Modified version that act as a function script and is to be used with Flow_Fit_ODE3_multi.m


% avaltouse is the a value to find R_b for  
% This data only works for Bo = 1
function Rb = Rb_Rc_fitting_func(avaltouse) 

    % Load data.  This data is for Bo=1, that's all we'll need for now.
    load('avec_Bo1_TWFamily.mat');
    load('Rminvec_Bo1_TWFamily.mat');
    % Find a-values closest to avaltouse and interpolate
    for jj=1:length(avec)-1
        if avec(jj+1)>avaltouse && avec(jj)<avaltouse
            break
        end
    end
    if jj<length(avec)-1
        Rmintouse=Rminvec(jj)+(Rminvec(jj+1)-Rminvec(jj))/(avec(jj+1)-avec(jj))*(avaltouse-avec(jj));
    else
        Rmintouse=NaN;
    end
    Rb = Rmintouse;
    
end
