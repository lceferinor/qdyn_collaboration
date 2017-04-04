%sp = FSEM3D_snapshot(0);
%figure,scatter(sp.X,sp.Z,[],sp.Tz)

function Steady_state(snap)
if nargin<2, cbar='off'; end
km=1000;
%snap=1015;
qdynfile=['fort.' num2str(snap)];

qdynsnap= Qdyn_read_ox_seq(qdynfile);

qdynfile0=['fort.20001'];
qdynsnap0= Qdyn_read_ox_seq(qdynfile0);


disp('Loading qdyn.in ...');
p=Qdyn_read_in('qdyn.in');
disp('Done');


qv=qdynsnap.v;
qth = qdynsnap.th;
qDc = p.DC;

%qdynpar = qv.*qth./qDc; % Gama, Gama=steady_state
%qdypar=qdynsnap.d;
%idst1=find(qdynpar<=0.99);
%qdynpar(idst1)=0;
%idst2=find(qdynpar>=1.001);
%qdynpar(idst2)=0;
Vdyninter=2*p.A.*qdynsnap.sigma./p.MU.*p.VS; % Gama
Vdynsteady=2*(p.A-p.B).*qdynsnap.sigma./p.MU.*p.VS;

%qdynpar=qv./Vdynsteady;
%qdynpar=p.DC;
qdynpar=qdynsnap.th;

fig=figure;scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdynpar)
year = 3600*24*365;

title([num2str((qdynsnap.t-qdynsnap0.t)) ' seconds']);

if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='northoutside';
%c.Label.String = 'Slip rate (m/s)';
%c.Label.String = 'Slip (m)';
c.Label.String = 'Gama';

end
caxis([0 max(qdynpar)])
%caxis([0 0.1]);
%caxis([0 0.01]);
%caxis([0 1]);

axis equal
xlim([min(qdynsnap.X/km) max(qdynsnap.X/km)]);
ylim([min(qdynsnap.Z/km) max(qdynsnap.Z/km)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
figname=['fig' num2str(snap)];
%figname=['fig_slip_qdyn' snap];
%print(fig,figname,'-djpeg')



end

