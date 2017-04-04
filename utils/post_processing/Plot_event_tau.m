%sp = FSEM3D_snapshot(0);
%figure,scatter(sp.X,sp.Z,[],sp.Tz)

function Plot_event_tau(snap,snap0)
if nargin<2, snap0=20001; end
km=1000;
%snap=1015;
qdynfile=['fort.' num2str(snap)];
qdynfile0=['fort.' num2str(snap0)];

qdynsnap= Qdyn_read_ox_seq(qdynfile);
qdynsnap0= Qdyn_read_ox_seq(qdynfile0);
%% dynamic parameters
p=Qdyn_read_in();
sigma = qdynsnap.sigma;
v=qdynsnap.v;
th=qdynsnap.th;
tau = sigma.*(p.MU_SS + p.A.*(log(v./p.V_SS)) + p.B.*(log(th.*p.V_SS./p.DC)));

qdypar=tau;
%qdypar=qdynsnap.d;

fig=figure;scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdypar/1e6);
year = 3600*24*365;

title([num2str((qdynsnap.t-qdynsnap0.t)/year) ' years']);

FigHandle = fig;
set(FigHandle, 'Position', [100, 100, 1049, 500]);
c=colorbar;
colormap jet
c.Location ='northoutside';
c.Label.String = 'Stress drop (MPa)';
caxis([0 max(qdypar/1e6)]);ylim([-25 0]);
xlim([min(qdynsnap.X/km) max(qdynsnap.X/km)]);
%ylim([min(p.Z/1e3) max(p.Z/1e3)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'color','w');
figname=['Stressdrop'];
print(FigHandle,figname,'-djpeg')




end

