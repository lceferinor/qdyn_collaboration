%sp = FSEM3D_snapshot(0);
%figure,scatter(sp.X,sp.Z,[],sp.Tz)

function [ti,ivmax]=Plot_event_tau_sliprate(snap,snap0)
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
ivmax=max(v);
th=qdynsnap.th;
tau = sigma.*(p.MU_SS + p.A.*(log(v./p.V_SS)) + p.B.*(log(th.*p.V_SS./p.DC)));
ti = qdynsnap.t(1);
qdypar=tau;
qdypar1=qdynsnap.v;

fig=figure;
subplot(2,1,1);scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdypar/1e6);
year = 3600*24*365;
title([num2str((ti)/year) ' years']);
FigHandle = fig;
set(FigHandle, 'Position', [100, 100, 1048, 1000]);
c=colorbar;
colormap jet
c.Location ='northoutside';
c.Label.String = 'Shear stress (MPa)';
xlim([min(qdynsnap.X/km) max(qdynsnap.X/km)]);
%caxis([0 max(qdypar/1e6)]);ylim([-25 0]);
caxis([0 50]);ylim([-25 0]);
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'color','w')
%figname=['Stressdrop'];
%print(FigHandle,figname,'-djpeg')

subplot(2,1,2);scatter(qdynsnap.X/km,qdynsnap.Z/km,[],qdypar1);

title([num2str((ti)/year) ' years']);
FigHandle = fig;
set(FigHandle, 'Position', [100, 100, 1100, 1200]);
c=colorbar;
colormap jet
c.Location ='northoutside';
c.Label.String = 'Slip velocity (m/s)';
caxis([0 0.1]);ylim([-25 0]);

xlim([min(qdynsnap.X/km) max(qdynsnap.X/km)]);
%ylim([min(p.Z/1e3) max(p.Z/1e3)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'color','w')


end

