%sp = FSEM3D_snapshot(0);
%figure,scatter(sp.X,sp.Z,[],sp.Tz)

function Steady_state(snap)
if nargin<2, cbar='off'; end
km=1000;
%snap=1015;
qdynfile=['fort.' num2str(snap)];

qdynsnap= Qdyn_read_ox_seq(qdynfile);
qdynsnap0= Qdyn_read_ox_seq(qdynfile0);
qdypar=qdynsnap.sigma;
%qdypar=qdynsnap.d;

fig=figure;scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdypar)
year = 3600*24*365;

title([num2str((qdynsnap.t-qdynsnap0.t)/year) ' years']);

if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='northoutside';
%c.Label.String = 'Slip rate (m/s)';
%c.Label.String = 'Slip (m)';
c.Label.String = 'Sigma (Pa)';

end
caxis([0 max(qdypar)])
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

