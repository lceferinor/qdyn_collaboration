
%%% A regime
fig=figure,
MoAdyn = [1.4564e20 1.3679e+20 1.5911e+20 1.5269e20 1.5531e+20];
RupAreaAdyn = [1.483e9 1.4749e+09 1.4988e+09 1.4928e9 1.9093e+09];

MoAqyn = [7.195e19];
RupAreaAqyn = [1.43e9];

loglog(MoAdyn,RupAreaAdyn,'or',MoAqyn,RupAreaAqyn,'ob')
hold on

%%% B regime (type II)
MoBdyn = [3.1102e19 3.22e19 2.2888e18  1.7748e18  3.2223e19];
RupAreaBdyn = [6.47e8 6.39343750e8 1.503125e8 140546875 639343750];

%3.5815e20 3.9451e20
%1.9354e9 1.9354e9
%3.593e20
%1.9354e9

MoBqyn = [1.693e19];
RupAreaBqyn = [4.993e8];

loglog(MoBdyn,RupAreaBdyn,'or',MoBqyn,RupAreaBqyn,'ob')
hold on

%%% C regime
MoCdyn = [1.0711e19 1.4262e20 1.3072e+18  6.2931e+17 5.7920e+18];
RupAreaCdyn = [3.38e8 1.9082e9 120703125 70890625 286125000] ; 

%2.6677e+17
%86578125
% 2.8799e+20
% 1.9354e+09
%3.0526e+20
%1.9354e+09

MoCqyn = [2.0475e18 ];
RupAreaCqyn = [1.355e8 ];

loglog(MoCdyn,RupAreaCdyn,'or',MoCqyn,RupAreaCqyn,'ob')
hold on

Moik=[1e17:1e17:5e20];

ik1=find(Moik<7.5e18);
Sik1=2.23*1e-15*(Moik(ik1)*1e7).^(2/3);
ik2=find(Moik>7.5e18);
Sik2=4.24*1e-11*(Moik(ik2)*1e7).^(1/2);

loglog(Moik(ik1),Sik1*1e6,'g',Moik(ik2),Sik2*1e6,'g')


%xlim([1e18 1e21])
%ylim([1e8 1e10])
legend('Fully-dynamic','Quasi-dynamic')
xlabel('Seismic moment (N.m)','FontSize',20)
ylabel('Rupture area (m^2)','FontSize',20)
set(gca,'FontSize',20)
figname=['Fig_Mo_irikura_Area'];
%figname=['fig_slip' timesnap];
print(fig,figname,'-djpeg')
