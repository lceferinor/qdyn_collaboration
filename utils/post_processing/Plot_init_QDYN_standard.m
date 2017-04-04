disp('Loading qdyn.in ...');
p=Qdyn_read_in_old('qdyn.in');
%p=Qdyn_read_in('qdyn.in');

disp('Done');

disp(['Loading QDYN output :  fort. ...']);
o = Qdyn_read_ox_seq('fort.1001');
disp('Done');
p.X=o.X;
p.Y=o.Y;
p.Z=o.Z;

qvar=p.DC;
%qvar=p.V_0;

fig=figure;scatter(p.X/1e3,p.Z/1e3,[],qvar)
FigHandle = fig;
set(FigHandle, 'Position', [100, 100, 1049, 500]);
c=colorbar;
c.Location ='northoutside';
c.Label.String = 'DC (m)';
%c.Label.String = 'Slip (m)';
%c.Label.String = 'Slip rate (m/s)';
%caxis([0 max(p.DC)])
caxis([0 0.05]);ylim([-25 0])
%caxis([0 max(qvar)]);

%axis equal
%axis tight


xlim([min(p.X/1e3) max(p.X/1e3)]);
%ylim([min(p.Z/1e3) max(p.Z/1e3)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
set(gcf,'color','w');
figname=['figDC'];
print(FigHandle,figname,'-djpeg')