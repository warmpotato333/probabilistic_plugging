% Reed Ogrosky

% Set a value to find R_b for.  
avaltouse=1.35;
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
% Make figure showing min of R
% If avaltouse>a_c (approx 1.377), then Rmintouse will be NaN
figure;
plot(avec,Rminvec,'k','Linewidth',1.25);
hold on;
plot(avec(592),Rminvec(592),'r.','Markersize',14);
plot(avaltouse,Rmintouse,'b.','Markersize',14);
text(avaltouse-0.135,Rmintouse,strcat(['if a=',num2str(avaltouse),', then R_b=',num2str(Rmintouse)]),'Color','blue');
text(1.32,0.67,'R_c=0.6691','Color','red');
xlabel('a');
ylabel('R_{min}');