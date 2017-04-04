
function SPECFEM_slip_inputpar(snaptime,cbar,L,corr_l,Dc_mean,Dc_min)

if nargin<2, cbar='off',L='64 km',corr_l='1 Km',Dc_mean='0.5', Dc_min='0.01'; end

stime = snaptime*4e-3;
snapspec = FSEM3D_snapshot(snaptime);
%fig=figure;scatter(snapspec.X,snapspec.Z,[],snapspec.Vx)
slip=sqrt(snapspec.Dx.^2 + snapspec.Dz.^2);
fig=figure;scatter(snapspec.X,snapspec.Z,[],slip)

%title(num2str(stime));




if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='eastoutside';
%c.Label.String = 'Slip rate (m/s)';
c.Label.String = 'Slip (m)';
end

timesnap=num2str(stime);
%% Calculating the magnitude
[Mw,RupArea]=Moment_ruptureArea(snaptime);
Mws=num2str(Mw);

%title(['Slip rate(m/s), time: ' timesnap ' s'],'FontSize',20)
title(['Mw=' Mws ', L=' L ' Lcorr=' corr_l ', Dcmean=' Dc_mean ', Dcsigma=' Dc_min],'FontSize',13)

caxis([0 max(slip)]);
%caxis([0 max(snapspec.Vx)])
%caxis([0 7])
%caxis([0 10])

axis equal
xlim([min(snapspec.X) max(snapspec.X)]);
ylim([min(snapspec.Z) max(snapspec.Z)])
xlabel('X (km)','FontSize',14)
ylabel('Z (km)','FontSize',14)
set(gca,'FontSize',14)
%figname=['fig' timesnap];
figname=['fig_magnitude' num2str(snaptime)];
print(fig,figname,'-dtiff')