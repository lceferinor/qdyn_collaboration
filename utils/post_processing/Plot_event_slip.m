%sp = FSEM3D_snapshot(0);
%figure,scatter(sp.X,sp.Z,[],sp.Tz)

function Plot_event_slip(snap,snap0,cbar)
if nargin<2, snap0=20001;cbar='off'; end
if nargin<3, cbar='off'; end
km=1000;
%snap=1015;
qdynfile=['fort.' num2str(snap)];
qdynfile0=['fort.' num2str(snap0)];

qdynsnap= Qdyn_read_ox_seq(qdynfile);
qdynsnap0= Qdyn_read_ox_seq(qdynfile0);
%qdypar=qdynsnap.v;
qdypar=qdynsnap.d-qdynsnap0.d;
fig=figure;scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdypar)

FigHandle = fig;
h=set(FigHandle, 'Position', [100, 100, 1049, 500]);
%title(num2str(qdynsnap.t-qdynsnap0.t));

if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='northoutside';
%c.Label.String = 'Slip rate (m/s)';
c.Label.String = 'Slip (m)';
end
caxis([0 max(qdypar)])
%caxis([0 0.1]);
%caxis([0 7]);

%axis equal
xlim([min(qdynsnap.X/km) max(qdynsnap.X/km)]);
ylim([min(qdynsnap.Z/km) max(qdynsnap.Z/km)])
xlabel('X (km)','FontSize',25)
ylabel('Z (km)','FontSize',25)
set(gca,'FontSize',25)
%figname=['fig' num2str(snap)];
set(gcf,'color','w')
figname=['fig_slip_qdyn' num2str(snap)];
print(fig,figname,'-djpeg')



end

