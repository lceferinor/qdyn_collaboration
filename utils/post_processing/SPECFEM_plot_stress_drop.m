
function SPECFEM_plot_stress_drop(snaptime,cbar)

if nargin<2, cbar='off'; end

stime = snaptime*4e-3;
snapspec = FSEM3D_snapshot(snaptime);
snapspec0 = FSEM3D_snapshot(100);

%fig=figure;scatter(snapspec.X,snapspec.Z,[],snapspec.Vx)
%slip=sqrt(snapspec.Dx.^2 + snapspec.Dz.^2);
Tstress=sqrt(snapspec.Tx.^2 + snapspec.Tz.^2);
Tstress0=sqrt(snapspec0.Tx.^2 + snapspec0.Tz.^2);
stressdrop=Tstress0-Tstress;
fig=figure;scatter(snapspec.X,snapspec.Z,[],stressdrop/1e6)

%title(num2str(stime));

if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='northoutside';
c.Label.String = 'Stress drop (Mpa)';
end
%c.Label.String = 'Slip (m)';
timesnap=num2str(stime);
%title(['Slip rate(m/s), time: ' timesnap ' s'],'FontSize',20)
%title(['time: ' timesnap ' s'],'FontSize',20)
%caxis([0 max(slip)]);
%caxis([0 max(stressdrop/1e6)])
caxis([0 10])
%caxis([0 0.001])
%caxis([0 1])


axis equal
xlim([min(snapspec.X) max(snapspec.X)]);
ylim([min(snapspec.Z) max(snapspec.Z)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
figname=['fig_stressdrop' num2str(snaptime)];
%figname=['fig_slip' timesnap];
print(fig,figname,'-djpeg')