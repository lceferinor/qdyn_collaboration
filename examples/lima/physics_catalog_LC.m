clc; clearvars; close all;

zoom_t=50;
scale_min=18;
scale_max=19.5;
scale_min_y=7.8;
scale_max_y=9.0;
% p = Qdyn_read_in();
% ot1 = Qdyn_read_ot('fort.18');

fname='LC_3DL620nx512W207.0552nw128';
%fname='LC_3DL800nx512W207.0552nw128_v4';
t0=0;			%starting time
t_inc=20000;		%time of earch run (in years)
event_all=[];
nseg=8;
Event_catalog = zeros(nseg,t_inc);


filename = [fname,'.mat' ];
% display(['Data   ' num2str(iirun) ' of ' num2str(n_run) ': ']);
display(['Processing ' filename ]);
load(filename);
close all
year=3600*24*365;
Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
t_ini = min(ot1.t/year);
t_fin = max(ot1.t/year);
%t_fin = 10^4;

v_th = Vdyn/10;            %threshold for seismic event    
slip_th = .1; % Slip threshold (for rupture boundaries)
%v_th = Vdyn/100;
% time correction
% ot1.t=ot1.t+(t0+(iirun-1)*t_inc)*year;
% ox1.t=ox1.t+(t0+(iirun-1)*t_inc)*year;
    
h2=figure('Units','inches', 'Position',[0 0 7 9]);

%subplot(6,1,[1:2])
%plot(ot1.t/year,ot1.p*p.W*p.MU);
%set(gca,'XTickLabel','');
%ylabel('Moment : (N.m)');
%xlim([t_ini t_fin]);


subplot(22,4,1:12);
iseis = 0;
count_s = 0; 
ii_start = 0;
n_out = p.N/(p.NXOUT*p.NWOUT);
d_temp=zeros(n_out,1);
a_temp=zeros(n_out,1);


for i= 1:1:numel(ot1.t)
    if iseis == 0 && ot1.v(i) >= v_th   % seismic event start
        iseis = 1;
        count_s = count_s + 1;
        ii_start = i;
        display(['Detected Seismic Event # ',num2str(count_s), ' at t = ', num2str(ot1.t(i)/year), ' year'])
    end
    if iseis == 1 && ot1.v(i) < v_th    % seismic event end
        iseis = 0;
        event(count_s,1) = ot1.t(ii_start);     %t start
        event(count_s,2) = ot1.t(i);     %t end
        event(count_s,3) = max(ot1.v(ii_start:1:i));        % vmax during the event
        event(count_s,4) = ot1.p(i)-ot1.p(ii_start); %potency
        event(count_s,5) = event(count_s,4)*p.MU*p.W;   % Moment
%        event(count_s,5) = event(count_s,4)*p.MU;      %3D
        event(count_s,6) = (2/3)*(log10(event(count_s,5))) - 6.0;	%M0
        display(['Duration: ' , num2str((event(count_s,2)-event(count_s,1))/60), 'm']);  
        id_ox1_start=find(ox1.t < ot1.t(ii_start),1,'last'); %Last index w/o EQ
        id_ox1_end=find(ox1.t > ot1.t(i),1,'first');%First index w/ EQ
        d_temp(:)=0;
        a_temp(:)=0;
        
        for is=1:1:n_out
            %id_start=find(ox1.v(is,id_ox1_start:1:id_ox1_end)> v_th,1,'first');
            %id_end=find(ox1.v(is,id_ox1_start:1:id_ox1_end)< v_th,1,'last');
            id_start=find(ox1.v(is,id_ox1_start:1:id_ox1_end)> v_th,1,'first');
            id_end=find(ox1.v(is,id_ox1_start:1:id_ox1_end)> v_th,1,'last');
            if isempty(id_start) || isempty(id_end)
                d_temp(is)=0;
                a_temp(is)=0;                
            else
                d_temp(is)=ox1.d(is,id_ox1_start-1+id_end)-ox1.d(is,id_ox1_start-1+id_start-1);
                a_temp(is)=1;
            end
        end
        if sum(a_temp) > 0
            %event(count_s,7) = sum(a_temp)/p.N*p.L*p.W;   % Rupture Area 
            event(count_s,7) = sum(a_temp)/n_out*p.L*p.W;   % Rupture Area
            %event(count_s,8) = sum(d_temp)/p.N*p.L*p.MU*p.W;   % Moment Seismic
            event(count_s,8) = sum(d_temp)/n_out*p.L*p.MU*p.W;   % Moment Seismic
            event(count_s,9) = (2/3)*(log10(event(count_s,8))) - 6.0;	%M0_seis
            vvmax_seis = max(ox1.v(:,id_ox1_start:1:id_ox1_end)'); 
            %idx_left=find(vvmax_seis>=v_th,1,'first');		%cell ID of rupture tip left
            %idx_right=find(vvmax_seis>=v_th,1,'last');          %cell ID of rupture tip left
            %event(count_s,10) = ox1.x(idx_left);		%rupture tip left location
            %event(count_s,11) = ox1.x(idx_right);         %rupture tip right location
            rupture_ind=find(vvmax_seis>=v_th);		%cell ID of rupture tip left
            event(count_s,10) = min(ox1.x(rupture_ind))-p.NXOUT/2*p.L/p.NX;		%rupture tip left location
            event(count_s,11) = max(ox1.x(rupture_ind))+p.NXOUT/2*p.L/p.NX;         %rupture tip right location     
            [vv_hypo idx_hypo] = max(ox1.v(:,id_ox1_start));	% find slip rate and cell ID of hypocenter
            event(count_s,12) = ox1.x(idx_hypo);         % Hypo location
            log_vel_maps(count_s,:) = log10(vvmax_seis/p.V_SS);
            slip_maps(count_s,:) = ox1.d(:,id_ox1_end)-ox1.d(:,id_ox1_start);	%slip maps of this event;
            % Rupture based on slip
            rupture_ind=find(slip_maps(count_s,:)>=slip_th);
            event(count_s,13) = min(ox1.x(rupture_ind))-p.NXOUT/2*p.L/p.NX;		%rupture tip left location
            event(count_s,14) = max(ox1.x(rupture_ind))+p.NXOUT/2*p.L/p.NX;         %rupture tip right location  
            event(count_s,15) = sum(ox1.dtau(a_temp == 1, id_ox1_end+1) - ...
                                ox1.dtau(a_temp == 1, max(1,id_ox1_start-1)))/1E6/sum(a_temp); % Average drop of dtau
            display(['Mw: ', num2str(event(count_s,9))]);
            display(['delta_tau: ', num2str(event(count_s,15)),'MPa']);
        else
           event(count_s,7:15) = 0;
           display(['Mw: ', num2str(event(count_s,9)),'. delta_tau: ', num2str(event(count_s,15)),'MPa']);
           slip_maps(count_s,:) = -1;
           log_vel_maps(count_s,:) = 0;
        end
    end
end


    
event(find(event(:,8)==0),:)=[];
count_s = length(event(:,1));
h = stem((event(:,1)+(event(:,2)-event(:,1))/2)/year,event(:,9),'--ok','MarkerFaceColor','r');
ylim([min(event(:,9))*0.9 max(event(:,9))*1.1]);
% ylim([4 8]);
xlim([t_ini t_fin]);
ylabel('Mw');
set(gca,'XtickLabel',[]);


subplot(22,4,13:24)
semilogy(ot1.t/year,ot1.v)
hold on 
semilogy(ot1.t/year,Vdyn*ones(size(ot1.t)),'r--');
ylabel('V: (m/s)');
legend('V_{max}','V_{dyn}');
xlim([t_ini t_fin]);
ylabel('V: (m/s)');
set(gca,'XtickLabel',[]);

subplot(22,4,25:48)
hold on;

for i = 1:count_s
    plot([event(i,1)/year,event(i,1)/year],...
        [event(i,10)/1000,event(i,11)/1000],'Color','k','LineWidth',2);
end

box on
x_min = min(p.X)/1000;
x_max = max(p.X)/1000;
xlim([t_ini t_fin]);
ylim([x_min x_max]);

ylabel('X (km)');
xlabel('Time: (years)');


%% Plot slip distributions

nx_out = p.NX/p.NXOUT;
nw_out = p.NW/p.NWOUT;


%subplot(21,4,[53,54,57,58,61,62,65,66,69,70,73,74]);

ind_sel = [];
for j=1:nw_out
    ind_sel = [ind_sel,(p.NWOUT*p.NX*(j-1) + (1:p.NXOUT:p.NX))];
end
dx = max(p.X)/p.NX/1000;
dy = max(p.Y)/p.NW/1000;
dz = abs(min(p.Z))/p.NW/1000;

y_min = min(p.Y)/10^3;
y_max = max(p.Y)/10^3;


slip_max = 5;

% Event 1
subplot(22,4,4+[53,57,61,65,69,73]);
event_id = 1;
scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],slip_maps(event_id,:),'s','filled');
caxis([0 slip_max]);
ylim([x_min x_max]);
xlim([y_min y_max+1]);
view(2);  
ylabel('X (km)');
xlabel('Y (km)');
box on
title(['Year: ', num2str(round(event(event_id,1)/year))],'FontWeight','normal');

% Event 2
subplot(22,4,4+1+[53,57,61,65,69,73]);
event_id = 2;
scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],slip_maps(event_id,:),'s','filled');
caxis([0 slip_max]);
ylim([x_min x_max]);
xlim([y_min y_max+1]);
view(2);  
set(gca,'YtickLabel',[]);
xlabel('Y (km)');
title(['Year: ', num2str(round(event(event_id,1)/year))],'FontWeight','normal');

% Event 3
subplot(22,4,4+2+[53,57,61,65,69,73]);
event_id = 3;
scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],slip_maps(event_id,:),'s','filled');
caxis([0 slip_max]);
ylim([x_min x_max]);
xlim([y_min y_max+1]);
view(2);  
set(gca,'YtickLabel',[]);
xlabel('Y (km)');
title(['Year: ', num2str(round(event(event_id,1)/year))],'FontWeight','normal');


% Event 4
subplot(22,4,4+3+[53,57,61,65,69,73]);
event_id = 4;
scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],slip_maps(event_id,:),'s','filled');
caxis([0 slip_max]);
ylim([x_min x_max]);
xlim([y_min y_max+1]);
view(2);  
set(gca,'YtickLabel',[]);
xlabel('Y (km)');
title(['Year: ', num2str(round(event(event_id,1)/year))],'FontWeight','normal');


subplot(22,4,4+[81:84]);
set(gca,'Visible','off')
caxis([0 slip_max]);
h = colorbar('North');
xlabel(h, 'Slip (m)');



print(h2,'-depsc2',[fname, 'physics_catalog.eps']);


%% Plot velocity distributions
% ln_v_min = -2;
% ln_v_max = 8;
% % Event 1
% subplot(21,4,[53,57,61,65,69,73]);
% event_id = 1;
% scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],log_vel_maps(event_id,:),'s','filled');
% caxis([ln_v_min ln_v_max]);
% ylim([x_min x_max]);
% xlim([y_min y_max+1]);
% view(2);  
% ylabel('X (km)');
% xlabel('Y (km)');
% 
% % Event 2
% subplot(21,4,1+[53,57,61,65,69,73]);
% event_id = 2;
% scatter3(p.Y(ind_sel)/1000+dy/2*p.NWOUT,p.X(ind_sel)/1000+dx/2*p.NXOUT,p.Z(ind_sel)/1000+dz/2*p.NWOUT,[],log_vel_maps(event_id,:),'s','filled');
% caxis([ln_v_min ln_v_max]);
% ylim([x_min x_max]);
% xlim([y_min y_max+1]);
% view(2);  
% set(gca,'YtickLabel',[]);
% xlabel('Y (km)');
% 
% 
% subplot(21,4,81:84);
% set(gca,'Visible','off')
% caxis([ln_v_min ln_v_max]);
% h = colorbar('North');
% xlabel(h, 'Slip (m)');
% 
% 



% 
% h5=figure(5);
% subplot(4,1,1:3)
% scatter(log10(event(:,8)),log10(event(:,7)))
% ylabel('Log10 Rupture Area : (m^2)');
% xlim([scale_min scale_max]);
% ylim([scale_min_y scale_max_y]);
% 
% subplot(4,1,4)
% hist(log10(event(:,8)),40)
% xlim([scale_min scale_max]);
% ylabel('Events');
% xlabel('Log10 Seismic Moment : (N.m)');
% print(h5,'-depsc2',[filename, '_M_A_scaling_zoom.eps']);
% 

save([filename, '_events.mat'],'event')



% Write file csv
filename = 'physics-based-catalog.csv';
csvwrite(filename, [event(:,1)/year, event(:,9:11)]);




% event_all = [event_all;event];
% 
% clear event ot1 ox1
% 
% 
% display(['All events events of ' fname ]);
% display(['Total time simulated ' num2str(t_inc)]);
% display(['Total number of events detected ' num2str(numel(event_all(:,1))) ]);
% 
% 
% 
% h4=figure(4);
% subplot(4,1,1:3)
% scatter(log10(event_all(:,8)),log10(event_all(:,7)))
% ylabel('Log10 Rupture Area : (m^2)');
% title(['Mid = ',num2str(median(log10(event_all(:,8)))), '   Mean = ', num2str(mean(log10(event_all(:,8))))]);
% 
% subplot(4,1,4)
% hist(log10(event_all(:,8)),100)
% ylabel('Events');
% xlabel('Log10 Seismic Moment : (N.m)');
% print(h4,'-depsc2',[fname, '_M_A_scaling_all_events.eps']);
% 
% 
% h5=figure(5);
% subplot(4,1,1:3)
% scatter(log10(event_all(:,8)),log10(event_all(:,7)))
% ylabel('Log10 Rupture Area : (m^2)');
% xlim([scale_min scale_max]);
% ylim([scale_min_y scale_max_y]);
% 
% subplot(4,1,4)
% hist(log10(event_all(:,8)),100)
% xlim([scale_min scale_max]);
% ylabel('Events');
% xlabel('Log10 Seismic Moment : (N.m)');
% print(h5,'-depsc2',[fname, '_M_A_scaling_zoom_all_events.eps']);
% save([fname, '_all_events.mat'],'event_all')
