clc;

fname='LC_3DL620nx512W207.0552nw128';
filename = [fname,'.mat' ];
load(filename);
close all

name='rupture';
mname=[name '.avi'];
year=3600*24*365;
vvmin=-2;  %min v/vpl
vvmax=8;   %max v/vpl
min_z =-60;
max_z =0;
t_const=1;  %time window for plotting v, if t_const=1, constant twin window, else, adjustable
t_win=10;  %time window for plotting v if t_const=1;
i_view=3;  %3 for 3d viw, 2 for 2d view
max_v_t=log10(max(ox1.v)./p.V_SS);
mean_v_t=log10(mean(ox1.v)./p.V_SS);
iplot=1;



vidObj = VideoWriter(mname);
vidObj.Quality = 100;
%vidObj.FrameRate = 10;
vidObj.FrameRate = 50;
open(vidObj);


for i=1:iplot:numel(ox1.t)
% for i=1:10
  disp(['Plotting ', num2str(i),'/',num2str(numel(ox1.t))]);
  h1=figure(1);
  subplot(6,6,[7:6:25]);
  temp_v=ox1.v(:,i);
  nx_out = p.NX/p.NXOUT;
  nw_out = p.NW/p.NWOUT;
  
  ind_sel = [];
  for j=1:nw_out
    ind_sel = [ind_sel,(p.NWOUT*p.NX*(j-1) + (1:p.NXOUT:p.NX))];
  end
   
  temp_v=reshape(temp_v,nx_out,nw_out);
  hold on;
  plot(log10(max(temp_v)./p.V_SS),p.Z(ind_sel(1:nx_out:end))/1000,'b');
  plot(log10(mean(temp_v)./p.V_SS),p.Z(ind_sel(1:nx_out:end))/1000,'g');
  legend('v_m_a_x','v_m_e_a_n');
  xlabel('log_1_0(V/V_p_l)');
  ylabel('Depth : km');
  xlim([vvmin vvmax]);
  ylim([min_z max_z]);
  hold off;
  subplot(6,6,[31:36]);
  hold on;

  plot(ox1.t/year,max_v_t,'b');
  plot(ox1.t/year,mean_v_t,'g');
  plot(ox1.t(i)/year,max_v_t(i),'bo');
  plot(ox1.t(i)/year,mean_v_t(i),'go');
  legend('v_m_a_x','v_m_e_a_n');
  ylabel('log_1_0(V/V_p_l)');
  xlabel('Time: years');
  ylim([vvmin vvmax]);
  if t_const == 1
      xlim([ox1.t(i)/year-t_win ox1.t(i)/year+t_win]);
  else
     if i<numel(ox1.t)
        xlim([ox1.t(i)-(ox1.t(i+1)-ox1.t(i))*5 ox1.t(i)+(ox1.t(i+1)-ox1.t(i))*5]/year);
     else
        xlim([ox1.t(i)-(ox1.t(i)-ox1.t(i-1))*5 ox1.t(i)+(ox1.t(i)-ox1.t(i-1))*5]/year); 
     end
  end
  hold off;
  subplot(6,6,[2:6 8:12 14:18 20:24 26:30]);
  dx = max(p.X)/p.NX/1000;
  dy = max(p.Y)/p.NW/1000;
  dz = abs(min(p.Z))/p.NW/1000;
  scatter3(p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],log10(ox1.v(:,i)/p.V_SS),'s','filled');
  caxis([-2 8]);
  colorbar('NorthOutside');
  if i_view == 3
      axis equal;
      view(i_view);
  else
      view(i_view);  
  end
  title(['time = ', num2str(ox1.t(i)/year,'%15.8f'),' year']);
  zlabel('Depth : km');
  xlabel('Location along-strike: km');
  ylabel('Location Y: km');
%  print(h1,'-djpeg','-r1000',[name, num2str(i), '.jpg']);
%  print(h1,'-dpdf',[name,'_m', num2str(i), '.pdf']);
%  mov(i)=getframe;
  writeVideo(vidObj, getframe(h1));
  clf(h1);
end
%movie2avi(mov,mname);
close(vidObj);