clearvars
clc

disp('Loading model output files')

%Load output files
load('ws_Pd_new_1_sensit_1.mat')
load('ws_global_new_1_sensit_1.mat')
load('ws_diagnostics_1_sensit_1.mat')
load('ws_Cd_new_1_sensit_1.mat')
load('frac_JLJR_diag.mat')

%% Load other things of the TM
[theta_P_P,theta_P_D, theta_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=...
    func_feeding_kernels_2(param);

load('matrix_nocorrection_01.mat'); %loads transport matrix
load('grid.mat');
load('config_data.mat');
load('boxes.mat')
load bins
load idx_wc

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
mon2= [31 28 31 30 31 30 31 31 30 31 30 31];
ticks = (1+cumsum(mon));
ticks(1) = 1;
labels = [{'J'} {'F'} {'M'} {'A'} {'M'} {'J'} {'J'} {'A'} {'S'} {'O'} {'N'} {'D'}];

xp = x;
xp(1)=0;
xp(end)=360;

Tyear=365*2*yrs+1-365*2:2:365*2*yrs;

zper=permute(z,[3 2 1]);
mprojection='robinson';%'winkel';%'bsam';
logf=1; %1 if in log scale 0 if normal
int_avg=0; %we do the average if 0. We integrate if 1
cmap=viridis(20);
depth_h=10; %depth idx
plot_fncs = function_plots_TMM();

deltaz=dznom;
deltaz_per=permute(deltaz,[3,2,1]);
deltaz_vec=deltaz(izBox); %depth vector in TM form
%% create/clean diagnostic

disp('Prparing 4D matrices of diagnostics and biomass')

% Calculate Global NPP and export
Export = s_flux_D_diag + s_flux_F_diag;
NPP_year_integrated_1=nansum(nansum(NPP2_diag.*volb)).*1e-18;
Export_year_integrated_1=sum(sum(Export(izBox==2,:).*volb(izBox==2)./deltaz_vec(izBox==2))).*1e-18;
Export_year_integrated_2=sum(sum(Export(izBox==7,:).*volb(izBox==7)./deltaz_vec(izBox==7))).*1e-18;

disp(['Global annual NPP: ', num2str(NPP_year_integrated_1),' PgC year^{-1}'])
disp(['Global annual export at 120m: ', num2str(Export_year_integrated_1),' PgC year^{-1}'])
disp(['Global annual export at 1080m: ', num2str(Export_year_integrated_2),' PgC year^{-1}'])

% Get diagnostics in 4D matrices,
% dimension are ([lon,lat,depth,time]),
% size(128, 64, 15, 365).
N_t=matrixToGrid(Nd, [], 'boxes.mat', 'grid.mat');
Export_t=matrixToGrid(Export, [], 'boxes.mat', 'grid.mat');
Export_D_t=matrixToGrid(s_flux_D_diag, [], 'boxes.mat', 'grid.mat');
Export_F_t=matrixToGrid(s_flux_F_diag, [], 'boxes.mat', 'grid.mat');
Export_F_small_t=matrixToGrid(s_flux_F_diag_small, [], 'boxes.mat', 'grid.mat');
NPP_t=matrixToGrid(NPP2_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
GPP_t=matrixToGrid(GPP2_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
pred_P_t=matrixToGrid(pred_P_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
pred_CP_t=matrixToGrid(pred_CP_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
pred_C_on_C_t=matrixToGrid(pred_C_on_C_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
morthtl_t=matrixToGrid(morthtl_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
reproduction_t=matrixToGrid(reproduction_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
pred_C_Dtot_t=matrixToGrid(pred_C_Dtot_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
FPP_t=matrixToGrid(FPP_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
DPP_t=matrixToGrid(DPP_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
remin_t=matrixToGrid(remin_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
NPPpico_t=matrixToGrid(NPP_pico_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
NPPnano_t=matrixToGrid(NPP_nano_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
NPPmicro_t=matrixToGrid(NPP_micro_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
frac_repr=squeeze(sum(reproduction_t(:,:,1:3,:),3)./sum(NPP_t(:,:,1:3,:),3));

%% Get Biomasses in 4D
% dimension are ([lon,lat,depth,time]),
% size(128, 64, 15, 365).

%biomass protists
Pdeach=zeros(nb,365,nbrP);
st1=1;
st2=nb;
for i=1:nbrP
   Pdeach(:,:,i)=Pd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Pd_t=matrixToGrid(sum(Pdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; 

%biomass copepods
Cdeach=zeros(nb,365,param.nbr_Ctot);
st1=1;
st2=nb;
for i=1:param.nbr_Ctot
   Cdeach(:,:,i)=Cd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Cd_t=matrixToGrid(sum(Cdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total copepod biomass


% biomass deadfalls
Ddeach=zeros(nb,365,param.nbr_D);
st1=1;
st2=nb;
for i=1:param.nbr_D
   Ddeach(:,:,i)=Dd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Dd_t=matrixToGrid(sum(Ddeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total deadfalls biomass

% biomass fecal pellets
Fdeach=zeros(nb,365,param.nbr_fp);
st1=1;
st2=nb;
for i=1:param.nbr_fp
   Fdeach(:,:,i)=Fd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Fd_t=matrixToGrid(sum(Fdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total fecal pellets biomass

%% Calculate size spectrum for each bin

disp('Calulating size spectrum for each bin (might take some time)')

zs=1;%3 %idx of the chosen surface layer we integrate biomass 

%Pdeach is protists biomass with dimensions [nb,365,number of protists], where nb is the TM
%vector size 52749 (the same goes for Cdeach for copepods)
%deltaz_per is the depth bin size in m (I wanted depth integrated biomass of the surface layers)

Pd_surf=zeros(nx,ny,365,nbrP);
for i=1:nbrP
    Pspec1=matrixToGrid(Pdeach(:,:,i), [], 'boxes.mat', 'grid.mat');
    Pd_surf(:,:,:,i)=squeeze(nansum(Pspec1(:,:,1:zs,:).*deltaz_per(:,:,1:zs),3));    
end

Cd_surf=zeros(nx,ny,365,length(param.Wvec));
for i=1:length(param.Wvec)
    Cspec1=matrixToGrid(Cdeach(:,:,i), [], 'boxes.mat', 'grid.mat');
    Cd_surf(:,:,:,i)=squeeze(nansum(Cspec1(:,:,1:zs,:).*deltaz_per(:,:,1:zs),3));    
end

%make community size spectrum
bin_spec=logspace(-7,4,12*2); %boundaries of the community size spectrum
bin_spec_mean=geomean([bin_spec(1:end-1);bin_spec(2:end)]); %mean size for each bin
Size_spec=zeros(nx,ny,365,length(bin_spec_mean)); %normalized size spectrum
Size_spec_sheldon=zeros(nx,ny,365,length(bin_spec_mean)); %size spectrum with total biomass in each bin
delta_sizebin=bin_spec(2:end)-bin_spec(1:end-1);
% get total biomass in each size bin
for i=1:length(bin_spec_mean)
    binP=find(param.V>bin_spec(i) & param.V<bin_spec(i+1));
    binC=find(param.Wvec>bin_spec(i) & param.Wvec<bin_spec(i+1));
    Size_spec_sheldon(:,:,:,i)=(nansum(Cd_surf(:,:,:,binC),4)+nansum(Pd_surf(:,:,:,binP),4));
    Size_spec(:,:,:,i)=(nansum(Cd_surf(:,:,:,binC),4)+nansum(Pd_surf(:,:,:,binP),4))./delta_sizebin(i);
end

clearvars Pspec1 Cspec1 binP binC

%get parameters of the size spectrum (coeficient and exponent)
%takes some time to run...
Size_spec_expo=NaN(nx,ny,365);
Size_spec_coeff=NaN(nx,ny,365);

bin_spec_log=log10(bin_spec_mean)';
for ii=1:nx
    for  jj=1:ny
        for dd=1:365
           if any(squeeze(Size_spec(ii,jj,dd,:)>0))
               
                pp = polyfit(bin_spec_log,log10(squeeze(Size_spec(ii,jj,dd,:))),1);           
                Size_spec_expo(ii,jj,dd)=pp(1);
                Size_spec_coeff(ii,jj,dd)=10.^pp(2);
           end
        end
    end
end

%%
% load tempday
%% Calculate average trophic level for all locations
% uncomment this section if TL needs to be recalculated, right now this
% ones is saved in advance as it takes quite long to run (loaded at the end of this section)

% TLmean=zeros(128,64);
% TLmean_day=zeros(128,64,365);
% 
% for dayd=1:365
% for xid=1:128
%     for yid=1:64
%         
%         if bathy(xid,yid,2)==1
% 
% frac_auto=frac_JLJR_diag(izBox==1 & ixBox==xid & iyBox==yid,:,dayd);
% Eflow_dep=zeros(54,54,2);
% biomvec_dep=zeros(2,54);
% for idep=1:2 %average over the 2 first depth ranges
% 
%     C=squeeze(Cdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
%     P=squeeze(Pdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
%     theta_P=theta_P_P;
% 
%     totsizevec=[param.V,param.Wvec'];
% 
%     % idx_prot=find(predsort_idx<=14);
% 
%     idx_parday=1:2:365*2;
%     T_inn=tempday(izBox==1 & ixBox==xid & iyBox==yid,idx_parday(dayd));%T(idx_1);
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%            % temperature effects
%            I_in=param.I'.*param.Q_I.^((T_inn-param.Tref)/10);
%            AN_in=param.alpha_N.*param.Q_N.^((T_inn-param.Tref)/10);
%            mu_max_in=param.mu_max.*param.Q_k.^((T_inn-18)/10);
%            AF_in=param.alpha_F.*1.5.^((T_inn-param.Tref)/10);
%            alpha_c=param.alpha_c'.*1.5.^((T_inn-param.Tref)/10);
%            VF_in=param.VF.*2.^((T_inn-18)/10);
% 
%     ingestion_each=zeros(length(param.Wvec),14+length(param.Wvec));
%     contribution_preys_to_cops=zeros(length(param.Wvec),14+length(param.Wvec));
%      for i=1:length(param.Wvec) %predators
%         for j=1:14+length(param.Wvec) %preys           
%             if j<=14 
%                 E_C_tot=(theta_cop_P(i,j).*P(j));
%             else
%                 E_C_tot=(theta_cop_cop(i,j-14).*C(j-14));       
%             end
%             F_lvl=(alpha_c(i).*E_C_tot)./(alpha_c(i).*E_C_tot+I_in(i));
%             ingestion_each(i,j)=I_in(i).*F_lvl.*C(i);             
%         end
%         ingestion_tot=sum(ingestion_each(i,:),2);
%         for j=1:14+length(param.Wvec) %preys 
%             contribution_preys_to_cops(i,j)=ingestion_each(i,j)./ingestion_tot;
%         end
%      end
% 
%     J_F=zeros(14,14);
%     contribution_preys_to_prots=zeros(14,14);
%      for i=1:14 %predators
%         for j=1:i-1%size(E_P_tot,1) %preys            
%             E_P_tot=P(j).*theta_P(i,j);
%             J_F(i,j)=(VF_in(i).*AF_in(i).*E_P_tot)./(AF_in(i).*E_P_tot+VF_in(i)).*P(i);
%         end
%         ingestion_tot=sum(J_F(i,:),2);
% 
%         for j=1:i-1%size(E_P_tot,1) %preys               
%             contribution_preys_to_prots(i,j)=J_F(i,j)./ingestion_tot;
%         end
% 
%      end
% 
%     J_F2=[J_F,zeros(14,length(param.Wvec))];
%     contribution_preys_to_prots_2=[contribution_preys_to_prots,zeros(14,length(param.Wvec))];
%     contributions_dep(:,:,idep)=[contribution_preys_to_prots_2;contribution_preys_to_cops];
%     Eflow_dep(:,:,idep)=[J_F2;ingestion_each];
% 
%     biomvec_dep(idep,:)=[P,C];
% 
% end
% 
% Eflow=mean(Eflow_dep,3); %flow of energy between compartments
% biomvec=mean(biomvec_dep,1); %mean biomass in the two depth layers
% contributions=mean(contributions_dep,3);
% 
% %to get the trophic level we need to organize all organisms relative to the
% %sze of the prey they eat, since they all have different predator-prey mass
% %ratios
% preysize=[log10(param.V)-log10(param.beta_P),...
%     log10(param.Wvec(param.ind_act)')-log10(param.beta_act),...
%     log10(param.Wvec(param.ind_pass)')-log10(param.beta_pass)];
% 
% [predsort, predsort_idx]=sort(preysize); %now predators are sorted according to the size of the prey they eat
% 
% %fraction that is preffered to eat by each organisms
% fracp=contributions'; % rows are preys and columns are predators
% 
% TL=ones(1,size(fracp,2));
% TLpart=zeros(1,size(fracp,2));
% 
% frac_auto_2=[frac_auto,zeros(1,length(param.Wvec))];
% 
% %rows are prey, columns are predator
% for i=predsort_idx%size(fracp,2) %predators
%     for j=1:i-1%size(E_P_tot,1) %preys
% 
%                 TLpart(j)=TL(j).*fracp(j,i).*(1-frac_auto_2(i));
% 
%     end
%     TL(i)=1+sum(TLpart);
% end
% 
% TLmean(xid,yid)=sum(TL(14+param.ind_act).*biomvec(14+param.ind_act))./(sum(biomvec(14+param.ind_act)));
% TLmean_day(xid,yid,dayd)=TLmean(xid,yid);
% 
%         end
%     end
% %     disp(xid)
% end
%     disp(dayd)
% end
% 
% TLmean_year=nanmean(TLmean_day,3);
% 
% save('TLmean_year.mat','TLmean_year');

load TLmean_year.mat


%% plots biomass spectrum global

x0=0;
y0=0;
width=18;
height=10;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

ax1=nexttile;
h=function_plot_map_contours(nanmean(nansum(Pd_t(:,:,1:3,:),3),4)./1000,0:9);
ylabel(h,'[gC m^{-2}]')
title('{\bfa.} Protists biomass','fontweight','normal','fontsize',10)
colormap(ax1,viridis(100));
caxis([0 9])

ax1=nexttile;
h=function_plot_map_contours(nanmean(nansum(Cd_t(:,:,1:3,:),3),4)./1000,0:9);
ylabel(h,'[gC m^{-2}]')
title('{\bfb.} Copepods biomass','fontweight','normal','fontsize',10)
colormap(ax1,viridis(100));
caxis([0 9])

ax2=nexttile(4);
h=function_plot_map_contours(mean(Size_spec_expo,3),-1.4:0.05:-0.6);%,0:0.2:2.8)
ylabel(h,'[-]')
caxis([-1.3 -0.7])
h.Ticks =-1.2:0.2:-0.6 ; 
title('{\bfd.} Exponent (\lambda)','fontweight','normal','fontsize',10)
colormap(ax2, brewermap(100,'BrBG'))

ax3=nexttile(3);
h=function_plot_map_contours(TLmean_year,2:0.25:4);
caxis([2 4])
ylabel(h,'[-]')
title('{\bfc.} Trophic level active copepods','fontweight','normal','fontsize',10)
colormap(ax3,viridis(100));

% print -depsc figures_paper/fig3_size_spec_global

%% biomass auto - mixo - hetero
%calculate the fraction of autotrophs, mixotrophs and heterotrophs
frac_JLJR_diag_per=permute(frac_JLJR_diag,[1,3,2]);

frac_auto=zeros(size(frac_JLJR_diag_per));
frac_mixo=zeros(size(frac_JLJR_diag_per));
frac_hetero=zeros(size(frac_JLJR_diag_per));

autolim=0.7; %if more than 0.7 of growth is driven by photo, then auto
heterolim=0.3;%if less than 0.3 of growth is driven by photo, then hetero
%in-betweens are considered mixos

frac_auto(frac_JLJR_diag_per>=autolim)=1;
frac_mixo(frac_JLJR_diag_per>heterolim & frac_JLJR_diag_per<autolim)=1;
frac_hetero(frac_JLJR_diag_per<=heterolim)=1;

biom_auto_t=matrixToGrid(sum(Pdeach.*frac_auto,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_mixo_t=matrixToGrid(sum(Pdeach.*frac_mixo,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_hetero_t=matrixToGrid(sum(Pdeach.*frac_hetero,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;

clearvars frac_JLJR_diag_per frac_auto frac_mixo frac_hetero

%% small cops biomass
% calculate small vs large copepod biomass
biom_cops_small_t=matrixToGrid(sum(Cdeach(:,:,param.Wvec<=1),3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_cops_large_t=Cd_t-biom_cops_small_t;

%% seasoanl
%plot seasonal dynamics

Export_F_large_t=Export_F_t-Export_F_small_t;

ordervec=[1:3:21,2:3:21,3:3:21,4:3:21];
nbrows=6;
nbcols=3;

st=1;
x0=0;
y0=0;
width=18;
height=15;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
left_color=[0 0 0];
right_color=[0.5 0.5 0.5];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
tiledlayout(7,3,'TileSpacing','Compact')%,'Padding','Compact');
% ha = tight_subplot(6,3,[.01 .05],[.1 .05],[.13 .12]);

for i=1:3
    if i==1
        x1=63;%62;
        y1=49;%51;
    elseif i==2
        x1=119;
        y1=52;
    elseif i==3
        x1=69;
        y1=40;
    end
    
cmap=brewermap(10,'Blues');
% subplot(nbrows,nbcols,ordervec(st))
nexttile(ordervec(st))
%   axes(ha(ordervec(st)))
plot(squeeze(nansum(Pd_t(x1,y1,1:3,:),3))./1000,'k','linewidth',1.2)
hold on
% plot(squeeze(nansum(Cd_t(x1,y1,1:3,:),3))./1000,'color',[0.5 0.5 0.5],'linewidth',1.2)

biom_prot_bar=[squeeze(nansum(biom_auto_t(x1,y1,1:3,:),3))./1000,...
    squeeze(nansum(biom_mixo_t(x1,y1,1:3,:),3))./1000,...
    squeeze(nansum(biom_hetero_t(x1,y1,1:3,:),3))./1000];
cmap=brewermap(3,'Blues');
h1=bar(biom_prot_bar,'stacked');
hold on
set(h1(1),'Facecolor',cmap(1,:)); 
set(h1(2),'Facecolor',cmap(2,:)); 
set(h1(3),'Facecolor',cmap(3,:)); 
set(h1,'FaceAlpha',0.6); 

xticklabels([])
if i==1
     title(' {\bfa} Protists','fontweight','normal')
elseif i==2
     title(' {\bfb}','fontweight','normal')
elseif i==3
     title(' {\bfc}','fontweight','normal')  
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend([h1(1),h1(2),h1(3)],'Auto.','Mixo.','Hetero.')
%     [lh,icons]=legend([h1(1),h1(2),h1(3)],'Auto.','Mixo.','Hetero.')
    legend boxoff
%     p1 = icons(1).Position;
% icons(1).Position = [0.4 p1(2) 0];
%     p1 = icons(2).Position;
% icons(2).Position = [0.4 p1(2) 0];
%     p1 = icons(3).Position;
% icons(3).Position = [0.4 p1(2) 0];
% hLegendPatch = findobj(icons, 'type', 'patch');
% set(hLegendPatch, 'XData', [.2, .2, .3, .3])

end
if i==1
%     ylabel({'Biomass','[mgC m^2]'})
    ylabel({'mgC m^{-2}'})
end


st=st+1;



nexttile(ordervec(st))
plot(squeeze(nansum(Cd_t(x1,y1,1:3,:),3))./1000,'color',[0.5 0.5 0.5],'linewidth',1.2)
hold on
biom_cop_bar=[squeeze(nansum(biom_cops_small_t(x1,y1,1:3,:),3))./1000,...
    squeeze(nansum(biom_cops_large_t(x1,y1,1:3,:),3))./1000];
cmap=brewermap(3,'Blues');
h1=bar(biom_cop_bar,'stacked');
hold on
set(h1(1),'Facecolor',[212, 158, 23]./255); 
set(h1(2),'Facecolor',[212, 158, 23]./(2.*255)); 
set(h1,'FaceAlpha',0.6); 
xticklabels([])
if i==1
     title(' {\bfd} Copepods','fontweight','normal')
elseif i==2
     title(' {\bfe}','fontweight','normal') 
elseif i==3
     title(' {\bff}','fontweight','normal')    
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend([h1(1),h1(2)],'Small','Large')
    legend boxoff
end
if i==1
    ylabel({'mgC m^{-2}'})
end


st=st+1;




nexttile(ordervec(st))
plot(squeeze(nansum(NPP_t(x1,y1,:,:),3)),'k','linewidth',1.2)
hold on
NPPsizes_bar=[squeeze(nansum(NPPpico_t(x1,y1,:,:),3)),...
    squeeze(nansum(NPPnano_t(x1,y1,:,:),3)),...
    squeeze(nansum(NPPmicro_t(x1,y1,:,:),3))];
cmap=brewermap(3,'Greens');
h1=bar(NPPsizes_bar,'stacked');
hold on
set(h1(1),'Facecolor',cmap(1,:)); 
set(h1(2),'Facecolor',cmap(2,:)); 
set(h1(3),'Facecolor',cmap(3,:)); 
set(h1,'FaceAlpha',0.6); 
xticklabels([])
if i==1
    title(' {\bfg} NPP','fontweight','normal')
elseif i==2
    title(' {\bfh}','fontweight','normal')    
elseif i==3
    title(' {\bfi}','fontweight','normal')        
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend(h1,'Pico','Nano','Micro')
    legend boxoff
end
if i==1
%     ylabel({'NPP','[mgC m^2 day^{-1}]'})
    ylabel({'mgC';'m^{-2} d^{-1}'})
end


st=st+1;
nexttile(ordervec(st))
exprt_bar_1=[squeeze(Export_D_t(x1,y1,2,:)),...
    squeeze(Export_F_small_t(x1,y1,2,:)),...
    squeeze(Export_F_large_t(x1,y1,2,:))];
cmap=brewermap(10,'Blues');
h1=bar(exprt_bar_1,'stacked');
hold on
set(h1(1),'Facecolor',cmap(end,:)); 
set(h1(2),'Facecolor',[212, 158, 23]./255); 
set(h1(3),'Facecolor',[212, 158, 23]./(2.*255)); 
set(h1,'FaceAlpha',0.6); 
h3=plot(squeeze(Export_t(x1,y1,2,:)),'k','linewidth',1.2);
xticklabels([])
if i==1
    title(' {\bfj} Export 120m','fontweight','normal')
elseif i==2
    title(' {\bfk}','fontweight','normal')    
elseif i==3
    title(' {\bfl}','fontweight','normal')        
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend([h1(1),h1(2),h1(3)],'Deadfall','FP_{small}','FP_{large}')
    legend boxoff
end
if i==1
    ylabel({'mgC';'m^{-2} d^{-1}'})
end

st=st+1;
nexttile(ordervec(st))
exprt_bar_2=[squeeze(Export_D_t(x1,y1,7,:)),...
    squeeze(Export_F_small_t(x1,y1,7,:)),...
    squeeze(Export_F_large_t(x1,y1,7,:))];
cmap=brewermap(10,'Blues');
hold on
h2=bar(exprt_bar_2,'stacked');
set(h2,'FaceAlpha',0.6); 
set(h2(1),'Facecolor',cmap(end,:)); 
set(h2(2),'Facecolor',[212, 158, 23]./255); 
set(h2(3),'Facecolor',[212, 158, 23]./(2.*255)); 
h4=plot(squeeze(Export_t(x1,y1,7,:)),'k:','linewidth',1.5);
xticklabels([])
if i==1
    title(' {\bfm} Export 1080m','fontweight','normal')
elseif i==2
    title(' {\bfn}','fontweight','normal')    
elseif i==3
    title(' {\bfo}','fontweight','normal')        
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend([h2(1),h2(2),h2(3)],'Deadfall','FP_{small}','FP_{large}')
    legend boxoff
end
if i==1
    ylabel({'mgC';'m^{-2} d^{-1}'})
end


    pe_ratio_day_100=zeros(365,1);
    pe_ratio_day_1000=zeros(365,1);
    for dmday=1:365
        if dmday<15
            nppavg=nanmean(cat(4,nansum(NPP_t(x1,y1,1:2,1:dmday),3),nansum(NPP_t(x1,y1,1:2,end-14+dmday:end),3)),4);
            pe_ratio_day_100(dmday)=squeeze(Export_t(x1,y1,2,dmday)./nppavg);
            pe_ratio_day_1000(dmday)=squeeze(Export_t(x1,y1,7,dmday)./nppavg);
        else
            pe_ratio_day_100(dmday)=squeeze(Export_t(x1,y1,2,dmday)./nanmean(nansum(NPP_t(x1,y1,1:2,dmday-14:dmday),3),4)); 
            pe_ratio_day_1000(dmday)=squeeze(Export_t(x1,y1,7,dmday)./nanmean(nansum(NPP_t(x1,y1,:,dmday-14:dmday),3),4));
        end
    end


st=st+1;
nexttile(ordervec(st))

hold on
yyaxis left
plot(pe_ratio_day_100,'-','color','k','linewidth',1.5)
hold on
if i==1
    ylabel({'[-]'})
end
if i==1
    title(' {\bfp} pe-ratio','fontweight','normal')
    ylim([0 0.5])
    yticks(0:0.2:0.5)
elseif i==2
    title(' {\bfq}','fontweight','normal')    
    ylim([0 0.5])
    yticks(0:0.2:0.5)
elseif i==3
    title(' {\bfr}','fontweight','normal')        
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';

yyaxis right
plot(pe_ratio_day_1000,'-','color',[0.5 0.5 0.5],'linewidth',1.5)
xticklabels([])
if i==3
yticks([4e-3 6e-3])
yticklabels({'0.004','0.006'})
end

set(gca,'box','off') 
xlim([0 365])
if i==3
    legend('120m','1000m')
    legend boxoff
end
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

corr(squeeze(Size_spec_expo(x1,y1,:)),squeeze(Export_t(x1,y1,2,:)./nansum(NPP_t(x1,y1,1:2,:),3)))
corr(squeeze(Size_spec_expo(x1,y1,:)),squeeze(Export_t(x1,y1,7,:)./nansum(NPP_t(x1,y1,:,:),3)))

corr(squeeze(Size_spec_expo(x1,y1,:)),pe_ratio_day_100)
corr(squeeze(Size_spec_expo(x1,y1,:)),pe_ratio_day_1000)


st=st+1;
nexttile(ordervec(st))
plot(squeeze(Size_spec_expo(x1,y1,:)),'k','linewidth',1.2)
xlabel('Days')
if i==1
    title(' {\bfs} Exponent size-spec.','fontweight','normal')
elseif i==2
    title(' {\bft}','fontweight','normal')    
elseif i==3
    title(' {\bfu}','fontweight','normal')        
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
xlim([0 365])
ylim([-1.4 -0.5])
if i==1
    ylabel({'[-]'})
end

st=st+1;

end

set(gcf,'color','w');

% print -depsc figures_paper/fig3_seasonal

%% size spectrum for a location only
idd=1;
x1=119;
y1=52;

deltaC=param.deltaC;
deltaC(end,:)=deltaC(end-1,:);
deltaC=deltaC(:);

c_ca=[29 41 81]./255;
c_cp=[87 160 211]./255;
c_p= [236 197 68]./255;

day_vec=[1 130 174 287];

st=1;
x0=0;
y0=0;
width=15;
height=8;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
for id=1:4
    
    idd=day_vec(id);
    
nexttile
h1=loglog(param.V,squeeze(Pd_surf(x1,y1,idd,:))./param.delta_V','color',c_p,'linewidth',1.5);
hold on
st1=1;
st2=param.nbr_stages;
for i=1:length(param.Wa)
    if i<=param.C_sp_act
        h2=plot(param.Wvec(st1:st2),squeeze(Cd_surf(x1,y1,idd,st1:st2))./deltaC(st1:st2),'color',c_ca,'linewidth',1.5) ;
    else
        h3=plot(param.Wvec(st1:st2),squeeze(Cd_surf(x1,y1,idd,st1:st2))./deltaC(st1:st2),'color',c_cp,'linewidth',1.5) ;
    end
    hold on
st1=st2+1;
st2=st2+param.nbr_stages;
end
h4=loglog(bin_spec_mean,100.*squeeze(Size_spec(x1,y1,idd,:)),'color',[0.5 0.5 0.5],'linewidth',1.2);
h5=loglog(bin_spec_mean,100.*Size_spec_coeff(x1,y1,idd).*bin_spec_mean.^(Size_spec_expo(x1,y1,idd)),'--','color',[0.5 0.5 0.5],'linewidth',1.2);
ylim([1e-5 1e13])
xlim([1e-7 1e4])

if id==1
    text(0.05,0.95,strcat(['{\bfa.} Day=',num2str(idd)]),'Units','normalized','fontsize',9)
    xticks([1e-6 1e-4 1e-2 1e0 1e2 1e4]);
    xticklabels({})
    yticks([1e-5 1e0 1e5 1e10]);
    yticklabels({'10^{-5}','10^{0}','10^{5}','10^{10}'})
    ylabel({'Biomass spectrum','[mgC m^{-2} \mugC^{-1}]'})
elseif id==2
    text(0.05,0.95,strcat(['{\bfb.} Day=',num2str(idd)]),'Units','normalized','fontsize',9)  
    xticks([1e-6 1e-4 1e-2 1e0 1e2 1e4]);
    xticklabels({})
    yticks([1e-5 1e0 1e5 1e10]);
    yticklabels({})
elseif id==3
    text(0.05,0.95,strcat(['{\bfc.} Day=',num2str(idd)]),'Units','normalized','fontsize',9) 
    xticks([1e-6 1e-4 1e-2 1e0 1e2 1e4]);
    xticklabels({'10^{-6}','10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'})
    yticks([1e-5 1e0 1e5 1e10]);
    yticklabels({'10^{-5}','10^{0}','10^{5}','10^{10}'})
    ylabel({'Biomass spectrum','[mgC m^{-2} \mugC^{-1}]'})
    xlabel('Body mass [\mugC]')
elseif id==4
    text(0.05,0.95,strcat(['{\bfd.} Day=',num2str(idd)]),'Units','normalized','fontsize',9)    
    xticks([1e-6 1e-4 1e-2 1e0 1e2 1e4]);
    xticklabels({'10^{-6}','10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'})
    yticks([1e-5 1e0 1e5 1e10]);
    yticklabels({})
    xlabel('Body mass [\mugC]')
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 

if id==2
    [hl,icons]=legend([h1,h2,h3,h4,h5],{'Protists','Active cops.','Passive cops.','Community\newlinespec.\times100','OLS fit'})
    legend boxoff
    for i_ic=1:5
    p1 = icons(i_ic).Position;
    icons(i_ic).Position = [0.22 p1(2) 0];
    icons(i_ic+5).XData = [0.05 0.2];
    icons(i_ic+10).XData = [0.05 0.2];
    end
end

end
set(gcf,'color','w');

% print -depsc figures_paper/fig3_spectrum.eps

%% generate food web plots

load tempday

nbc=3;

st=1;
x0=0;
y0=0;
width=18;
height=6;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

for iloc=1:3
    
%     if iloc==1
%         xid=62;%117;
%         yid=51;
%         dayd=250;
%         
%     elseif iloc==2
%         
%         xid=62;%117;
%         yid=51;
%         dayd=250;
%         
%     elseif iloc==3
%         
%         xid=62;%117;
%         yid=51;
%         dayd=250;
%         
%     elseif iloc==4
%         
%         xid=62;%117;
%         yid=51;
%         dayd=250;
%         
%     end
    
    if iloc==1
        xid=63;%117;
        yid=40;
        dayd=200;%300;%135;
        
    elseif iloc==2
        
        xid=63;%117;
        yid=44;
        dayd=200;%300;%135;
        
    elseif iloc==3
        
        xid=63;%117;
        yid=51;
        dayd=200;%300;%135;
        
    end
    
% % xid=76;%65
% % yid=41;%49
% % % 
% % % xid=65;
% % % yid=49;
% % % dayd=130;
% % % % % % 
% % xid=119; %NA
% % yid=52;
% % dayd=350;
% % % 
% % xid=62; %NPA
% % yid=51;
% % dayd=300;
% % % 
% % xid=89;
% % yid=31;
% % dayd=300;



% dayd=287;
% xid=119;
% yid=52;

c_ca=[29 41 81]./255;
c_cp=[87 160 211]./255;
c_p= [236 197 68]./255;

day_vec=[1 130 174 287];


frac_auto=frac_JLJR_diag(izBox==1 & ixBox==xid & iyBox==yid,:,dayd);
Eflow_dep=zeros(54,54,2);
biomvec_dep=zeros(2,54);
for idep=1:2

C=squeeze(Cdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
P=squeeze(Pdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
theta_P=theta_P_P;

totsizevec=[param.V,param.Wvec'];

% idx_prot=find(predsort_idx<=14);

idx_parday=1:2:365*2;
T_inn=tempday(izBox==1 & ixBox==xid & iyBox==yid,idx_parday(dayd));%T(idx_1);

        %%%%%%%%%%%%%%%%%%%%%%%%%
       % temperature effects
       I_in=param.I'.*param.Q_I.^((T_inn-param.Tref)/10);
       AN_in=param.alpha_N.*param.Q_N.^((T_inn-param.Tref)/10);
       mu_max_in=param.mu_max.*param.Q_k.^((T_inn-18)/10);
       AF_in=param.alpha_F.*1.5.^((T_inn-param.Tref)/10);
       alpha_c=param.alpha_c'.*1.5.^((T_inn-param.Tref)/10);
       VF_in=param.VF.*2.^((T_inn-18)/10);
 
nu=zeros(length(param.Wvec),14+length(param.Wvec));
 for i=1:length(param.Wvec) %predators
    for j=1:14+length(param.Wvec) %preys           
        if j<=14 
            E_C_tot=(theta_cop_P(i,j).*P(j));
        else
            E_C_tot=(theta_cop_cop(i,j-14).*C(j-14));       
        end
        F_lvl=(alpha_c(i).*E_C_tot)./(alpha_c(i).*E_C_tot+I_in(i));
        nu(i,j)=I_in(i).*F_lvl.*C(i); 
    end
 end
 
J_F=zeros(14,14);
 for i=1:14 %predators
    for j=1:i-1%size(E_P_tot,1) %preys            
        E_P_tot=P(j).*theta_P(i,j);
        J_F(i,j)=(VF_in(i).*AF_in(i).*E_P_tot)./(AF_in(i).*E_P_tot+VF_in(i)).*P(i);
    end
 end

J_F2=[J_F,zeros(14,length(param.Wvec))];
Eflow_dep(:,:,idep)=[J_F2;nu];

biomvec_dep(idep,:)=[P,C];

end

Eflow=mean(Eflow_dep,3);
biomvec=mean(biomvec_dep,1);

preysize=[log10(param.V)-log10(param.beta_P),...
    log10(param.Wvec(param.ind_act)')-log10(param.beta_act),...
    log10(param.Wvec(param.ind_pass)')-log10(param.beta_pass)];

[predsort, predsort_idx]=sort(preysize);

frac_c=[(theta_cop_P.*P),(theta_cop_cop.*C)];
frac_p_to_c=zeros(size((theta_cop_P.*P(1,:))'));
frac_p=[P.*theta_P,frac_p_to_c];
fracp=[frac_p;frac_c]';
fracp=fracp(predsort_idx,predsort_idx);

TL=ones(1,size(fracp,2));
TLpart=zeros(1,size(fracp,2));

%rows are prey, columns are predator
for i=1:size(fracp,2) %predators
    for j=1:i-1%size(E_P_tot,1) %preys
        if i<15
%             if frac_auto(i)<0.90
                TLpart(j)=TL(j).*(fracp(j,i)./sum(fracp(:,i),1)).*(1-frac_auto(i));
%             else
%                 TLpart(j)=0;
%             end        
        else
            TLpart(j)=TL(j).*(fracp(j,i)./sum(fracp(:,i),1));        
        end
    end
    TL(i)=1+sum(TLpart);
end
% end

idx_prots=find(predsort_idx<=14);
idx_act=find(predsort_idx>=min(param.ind_act)+14 & predsort_idx<=max(param.ind_act)+14);
idx_pass=find(predsort_idx>max(param.ind_act)+14);

% subplot(1,nbc,iloc)
nexttile
Eflow_sort=Eflow(predsort_idx,predsort_idx);
clight=[237, 232, 218]./255;
for i=2:size(fracp,2) %predators
%     if max(theta_P(i,:)) > 1e-4
        for j=1:i-1%size(E_P_tot,1) %preys 
%             interc=(fracp(j,i)./sum(fracp(:,i),1));
            interc=(Eflow_sort(i,j))./sum(sum(Eflow));
%             if interc>1e-3 & biomvec(predsort_idx(i))./(sum(biomvec))>1e-2
            if interc>1e-3%>0.00001
                plot([totsizevec(predsort_idx(i)),totsizevec(predsort_idx(j))],[TL(i),TL(j)],'color',clight,'linewidth',interc.*30)
                hold on
            end
        end
end

bnorm=5000;
scatter(totsizevec(predsort_idx),TL,biomvec(predsort_idx)./bnorm.*1000,...
    'filled','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');%,'markerfacealpha',0.5)
% scatter(totsizevec(predsort_idx),TL,20.*(log10(biomvec(predsort_idx))./log10((sum(biomvec)))+abs(min(log10(biomvec(predsort_idx))./log10((sum(biomvec)))))+1e-3),...
%     'filled','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k')%,'markerfacealpha',0.5)
h1=scatter(param.V,TL(idx_prots),biomvec(predsort_idx(idx_prots))./bnorm.*1000,...
    'filled','markerfacecolor',c_p,'markeredgecolor','k');
h2=scatter(param.Wvec(predsort_idx(idx_act)-14),TL(idx_act),biomvec(predsort_idx(idx_act))./bnorm.*1000,...
    'filled','markerfacecolor',c_ca,'markeredgecolor','k');
h3=scatter(param.Wvec(predsort_idx(idx_pass)-14),TL(idx_pass),biomvec(predsort_idx(idx_pass))./bnorm.*1000,...
    'filled','markerfacecolor',c_cp,'markeredgecolor','k');
set(gca,'xscale','log')
ylim([0.5 6])
xlim([1e-8 1e5])
xlabel('Body mass [\mugC]')
xticks(logspace(-7,5,5));
set(gcf,'color','w');
if iloc==1
    text(0.05,0.98,strcat(['{\bfa.} Lat.=',num2str(y(yid),'%.0f')]),'Units','normalized','fontsize',10)
%     title(strcat(['   {\bfa} Lat.=',num2str(y(yid))]),'fontweight','normal')
ylabel('Trophic level')
elseif iloc==2
%     title(strcat(['   {\bfb} Lat.=',num2str(y(yid))]),'fontweight','normal') 
    text(0.05,0.98,strcat(['{\bfb.} Lat.=',num2str(y(yid),'%.0f')]),'Units','normalized','fontsize',10)
elseif iloc==3
%     title(strcat(['   {\bfc} Lat.=',num2str(y(yid))]),'fontweight','normal') 
    text(0.05,0.98,strcat(['{\bfc.} Lat.=',num2str(y(yid),'%.0f')]),'Units','normalized','fontsize',10)
end
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 

TLmean=sum(TL.*biomvec(predsort_idx))./(sum(biomvec));
TLmean_lvl1=sum((totsizevec(predsort_idx(TL<=2))).*biomvec(predsort_idx(TL<=2)))./(sum(biomvec(TL<=2)));

end

legend([h1,h2,h3],'Protists','Active copepods','Passive copepods')
legend boxoff

% print -depsc figures_paper/fig3_foodweb

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NPP field data compiled in Saba et al 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NPP_int=squeeze(nansum(NPP_t,3));

NPP_int(NPP_int<1)=NaN;

%Load dataset from Saba et al 2011
dat_2=importdata('data_PP_field/bg-8-489-2011-supplement/File_S1_PPARR4_Database.xlsx'); %Saba dataset

lat_1=dat_2.data(:,4);
lon_1=dat_2.data(:,5);
lon360 = wrapTo180(lon_1);
% lon_2=mod(lon_2(:),360);
day_1=dat_2.data(:,7);
PP_1=dat_2.data(:,13); %mgC m^-2 d^-1
%remove data from mediterranean
idx_med=find(lon360>0 & lon360<42 & lat_1>30 & lat_1<46);

lat_1(idx_med)=NaN;
lon_1(idx_med)=NaN;
day_1(idx_med)=NaN;
PP_1(idx_med)=NaN;

NPP_tmm=NaN(length(x),length(y),365);
NPP_tmm_std=NaN(length(x),length(y),365);
for i=1:length(x)-1
    for j=1:length(y)-1
        int_x=[x(i), x(i+1)];
        int_y=[y(j), y(j+1)];  
        idx_lonlat=find(lon_1>=int_x(1) & lon_1<int_x(2) & lat_1>=int_y(1) & lat_1<int_y(2));
        if ~isempty(idx_lonlat)
            for k=1:365
                idx_lonlat=find(lon_1>=int_x(1) & lon_1<int_x(2) & lat_1>=int_y(1) & lat_1<int_y(2) & day_1==k);
                if ~isempty(idx_lonlat)
                    NPP_tmm(i,j,k)=nanmean(PP_1(idx_lonlat));
                    NPP_tmm_std(i,j,k)=std(PP_1(idx_lonlat));
                end
            end
        end       
    end
end

lat_3D=repmat(y',[nx,1,365]);
lon_3D=repmat(wrapTo180(x),[1,ny,365]);

% figure
% subplot(2,2,1)
% scatter(NPP_tmm(NPP_tmm_std<0.1.*NPP_tmm),NPP_int(NPP_tmm_std<0.1.*NPP_tmm),10,lat_3D(NPP_tmm_std<0.1.*NPP_tmm));
% hold on
% plot(linspace(0,1500),linspace(0,1500),'k--')
% colorbar
% xlabel('Data [mgC m^{-2} day^{-1}]')
% ylabel('Model [mgC m^{-2} day^{-1}]')
% title('NPP')
% 
% subplot(2,2,2)
% histogram(log10(NPP_tmm(NPP_tmm_std<0.1.*NPP_tmm)));
% hold on
% histogram(log10(NPP_int(~isnan(NPP_tmm) & NPP_tmm_std<0.1.*NPP_tmm)));
% legend('Data','Model')
% xlabel('Log_{10} NPP [mgC m^{-2} day^{-1}]')
% ylabel('[]')
% title('Log_{10} NPP')
% 
% subplot(2,1,2)
% scatter(lon_3D(NPP_tmm_std<0.1.*NPP_tmm),lat_3D(NPP_tmm_std<0.1.*NPP_tmm),8,'g','filled');
% hold on
% geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
% title('NPP data locations')
% set(gcf,'color','w');

idx_SO=find(lat_3D<-50);
idx_HOT=find(lon_3D<-150 & lat_3D>0);
idx_BATS=find(lon_3D>-67 & lon_3D<-65 & lat_3D>29 & lat_3D<30);
idx_NA=find(lon_3D>-21 & lon_3D<0);
idx_arab=find(lon_3D>50 & lon_3D<67);


% NPP_int_log=log10(NPP_int(:));
% NPP_tmm_log=log10(NPP_tmm(:));
% NPP_int_log(isinf(NPP_int_log))=NaN;
% NPP_tmm_log(isinf(NPP_int_log))=NaN;
% NPP_int_log(isinf(NPP_tmm_log))=NaN;
% NPP_tmm_log(isinf(NPP_tmm_log))=NaN;
% 
% NPPso_nan=NPP_tmm(idx_SO);
% NPPso_nan=NPPso_nan(~isnan(NPPso_nan));
% NPPhot_nan=NPP_tmm(idx_HOT);
% NPPhot_nan=NPPhot_nan(~isnan(NPPhot_nan));
% NPPbats_nan=NPP_tmm(idx_BATS);
% NPPbats_nan=NPPbats_nan(~isnan(NPPbats_nan));
% NPPna_nan=NPP_tmm(idx_NA);
% NPPna_nan=NPPna_nan(~isnan(NPPna_nan));
% NPPar_nan=NPP_tmm(idx_arab);
% NPPar_nan=NPPar_nan(~isnan(NPPar_nan));
% 
% RMSD1=(1./length(NPPso_nan).*(nansum((log10(NPP_int(idx_SO))-log10(NPP_tmm(idx_SO))).^2))).^0.5;
% RMSD2=(1./length(NPPhot_nan).*(nansum((log10(NPP_int(idx_HOT))-log10(NPP_tmm(idx_HOT))).^2))).^0.5;
% RMSD3=(1./length(NPPbats_nan).*(nansum((log10(NPP_int(idx_BATS))-log10(NPP_tmm(idx_BATS))).^2))).^0.5;
% RMSD4=(1./length(NPPna_nan).*(nansum((log10(NPP_int(idx_NA))-log10(NPP_tmm(idx_NA))).^2))).^0.5;
% RMSD5=(1./length(NPPar_nan).*(nansum((log10(NPP_int(idx_arab))-log10(NPP_tmm(idx_arab))).^2))).^0.5;
% RMSD=[RMSD1;RMSD2;RMSD3;RMSD4;RMSD5];


% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CARBON EXPORT - DATA FROM HENSON 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% dat=importdata('data_carbon_export/gbc20895-sup-0002-2018gb006158-data set s1.xlsx');
% 
% lat_data_flux=dat.data(:,1);
% lon_data_flux=dat.data(:,2);
% day_data_flux=dat.data(:,3);
% flux_data=dat.data(:,4);
% PP_data_flux=dat.data(:,5);
% eratio_data=dat.data(:,6);
% depth_data=dat.data(:,7);
% macrozoo_data_flux=dat.data(:,15);
% mesozoo_data_flux=dat.data(:,13);
% 
% lon2_flux=mod(lon_data_flux(:),360);
% 
% zface=cumsum(deltaz);
% 
% mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
% monticks = [(1+cumsum(mon)) 365];
% 
% NPP_t_100=squeeze(nansum(NPP_t(:,:,1:2,:),3));
% 
% flux_model=zeros(length(lon_data_flux),1);
% eff_model=zeros(length(lon_data_flux),1);
% PP_model_flux=zeros(length(lon_data_flux),1);
% for i=1:length(lon_data_flux)
%     
%     if isnan(lat_data_flux(i)) || isnan(lat_data_flux(i))
%         xid=find(x>lon2_flux(i) & x<lon2_flux(i)+2.9);
%         yid=find(y>lat_data_flux(i) & y<lat_data_flux(i)+2.9);
%         flux_model(i)=NaN;
%         eff_model(i)=NaN;
%         PP_model_flux(i)=NaN;
%     else
%         xid=find(x>lon2_flux(i) & x<lon2_flux(i)+2.9);
%         if lon2_flux(i)>max(x)
%             xid=128;
%         end
%         yid=find(y>lat_data_flux(i) & y<lat_data_flux(i)+2.9);
%         if depth_data(i)<130
%             flux_model(i)=nanmean(Export_t(xid(1),yid(1),2,day_data_flux(i)),4);
%         elseif depth_data(i)>=130
%             flux_model(i)=nanmean(Export_t(xid(1),yid(1),3,day_data_flux(i)),4);
%         end
%         eff_model(i)=flux_model(i)./nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
%         PP_model_flux(i)=nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
%     end
% 
%     
% end
% 
% 
% idx_data=find(~isnan(eff_model) & ~isinf(eff_model) & ~isnan(eratio_data) & PP_model_flux>=0 & flux_data>0);
% 
% flux_model2=flux_model(idx_data);
% eff_model2=eff_model(idx_data);
% PP_model_flux2=PP_model_flux(idx_data);
% flux_data2=flux_data(idx_data);
% eratio2=eratio_data(idx_data);
% lat_data2=lat_data_flux(idx_data);
% lon_data2=lon_data_flux(idx_data);
% PP_data_flux2=PP_data_flux(idx_data);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export data from Lutz et al 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%model interpolated


filename='data_carbon_export/Lutz_2007/data_carbon_export_lutz_2007.xlsx';
lutzdata=importdata(filename);

lat_flux_lut=lutzdata.data(:,1);
lon_flux_lut_180=lutzdata.data(:,2);
lon_flux_lut=wrapTo360(lutzdata.data(:,2));
depth_flux_lut=lutzdata.data(:,3);
day_start_lut=lutzdata.data(:,6);
mon_start_lut=lutzdata.data(:,7);
year_start_lut=lutzdata.data(:,8);
day_end_lut=lutzdata.data(:,10);
mon_end_lut=lutzdata.data(:,11);
year_end_lut=lutzdata.data(:,12);
flux_data_lut=lutzdata.data(:,13);


z_interp=1:6200;

Flux_model_lut=NaN(length(lat_flux_lut),1);
Flux_model_lut_depth=NaN(length(lat_flux_lut),15);
Flux_model_lut_depth_interp=NaN(length(lat_flux_lut),length(z_interp));
zface=cumsum(deltaz);
for i=1:length(lat_flux_lut)
    idx_x=find(x>=lon_flux_lut(i)-3 & x<lon_flux_lut(i)+3);
    idx_y=find(y>=lat_flux_lut(i)-2.8 & y<lat_flux_lut(i)+2.8);
    idx_day_start=30.*mon_start_lut(i)-30+day_start_lut(i);
    idx_day_end=30.*mon_end_lut(i)-30+day_end_lut(i);
    if depth_flux_lut(i)<max(zface)
        idx_z=find(zface>depth_flux_lut(i));
    else
        idx_z=15;
    end
    if ~isempty(idx_x) && ~isnan(depth_flux_lut(i))
        if year_end_lut(i)-year_start_lut(i) <= 1 && mon_end_lut(i)-mon_start_lut(i) >=-1 && mon_end_lut(i)-mon_start_lut(i) <=1
          F_wc=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,:),4));
%           F_interp = interp1(zface,F_wc,z_interp);
% %           interp1(zface,F_wc,z_interp,'linear','extrap')
%           idx_z=find(z_interp>depth_flux_lut(i));          
%           Flux_model_lut(i)=F_interp(idx_z(1));
%               if isnan(F_interp(idx_z(1)))                  
%                   Flux_model_lut(i)=interp1(zface(~isnan(F_wc)),F_wc(~isnan(F_wc)),depth_flux_lut(i),'linear','extrap');
%               end
          Flux_model_lut_depth(i,:)=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,:),4))';
        elseif year_end_lut(i)-year_start_lut(i) <= 1 
          F_wc=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,idx_day_start:idx_day_end),4));
%           F_interp = interp1(zface,F_wc,z_interp);  
%           idx_z=find(z_interp>depth_flux_lut(i));
%           Flux_model_lut(i)=F_interp(idx_z(1));  
%               if isnan(F_interp(idx_z(1)))                  
%                   Flux_model_lut(i)=interp1(zface(~isnan(F_wc)),F_wc(~isnan(F_wc)),depth_flux_lut(i),'linear','extrap');
%               end
          Flux_model_lut_depth(i,:)=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,idx_day_start:idx_day_end),4))';
        elseif year_end_lut(i)-year_start_lut(i) ~=0 && mon_end_lut(i)-mon_start_lut(i) >=-1 && mon_end_lut(i)-mon_start_lut(i) <=1
          F_wc=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,:),4));
         Flux_model_lut_depth(i,:)=squeeze(nanmean(Export_t(idx_x(1),idx_y(1),:,:),4))'; 
        end
        
        if any(isfinite(F_wc))
        F_interp = interp1(zface,F_wc,z_interp);  
        idx_z=find(z_interp>depth_flux_lut(i));
        Flux_model_lut(i)=F_interp(idx_z(1));
        Flux_model_lut_depth_interp(i,:)=F_interp;
        
%        if depth_flux_lut(i)>5500
%            alk
%        end
        
              if isnan(F_interp(idx_z(1)))                  
                  Flux_model_lut(i)=10.^interp1(zface(~isnan(F_wc)),log10(F_wc(~isnan(F_wc))),depth_flux_lut(i),'linear','extrap');
                  F_interp(idx_z(1))=Flux_model_lut(i);
                  Flux_model_lut_depth_interp(i,:)=10.^interp1(z_interp(~isnan(F_interp)),log10(F_interp(~isnan(F_interp))),z_interp);
                  
%                   Flux_model_lut_depth_interp(i,:)=F_interp;
                  if Flux_model_lut(i)<0
                      Flux_model_lut(i)=NaN;
                  end
              end
        end 
        
        
        
        
        size(Flux_model_lut_depth)
    end
end



%% get Ocean basins for Lutz
%1- polar north
% 2:21 - Atlantic ocean
% 22:30 - Indian ocean
% 31:50 - Pacific ocean
% 50:54 - Southern ocean

s = shaperead('../3Doutputs/Longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp');

idx_polar=[32,31,3,2,1]; 
idx_temp=[30,29,20,18,17,16,13,5,4]; 
idx_trop=[26,25,14,9,8,27]; 
idx_oligo=[24,23,21,15,10,7,6];

%remove coastal regions
idx_open=zeros(1,54);
for i=1:54
   alk=split(s(i).ProvDescr);
   if ~contains(alk{1},'Coastal')
      idx_open(i)=i;
   end   
end
idx_open=idx_open(idx_open~=0);

% figure
% mapshow(s(idx_open(idx_trop)))
% hold on
% mapshow(s(idx_open(idx_oligo)))
% mapshow(s(idx_open(idx_temp)))
% mapshow(s(idx_open(idx_polar)))

idx_lutz_oligo=zeros(size(Flux_model_lut_depth));
for i=idx_open(idx_oligo)
    in = inpolygon(lon_flux_lut_180,lat_flux_lut,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_lutz_oligo(in==1)=1;
end
idx_lutz_oligo=logical(idx_lutz_oligo);

idx_lutz_polar=zeros(size(Flux_model_lut_depth));
for i=idx_open(idx_polar)
    in = inpolygon(lon_flux_lut_180,lat_flux_lut,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_lutz_polar(in==1)=1;
end
idx_lutz_polar=logical(idx_lutz_polar);

idx_lutz_temp=zeros(size(Flux_model_lut_depth));
for i=idx_open(idx_temp)
   in = inpolygon(lon_flux_lut_180,lat_flux_lut,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_lutz_temp(in==1)=1;
end
idx_lutz_temp=logical(idx_lutz_temp);

idx_lutz_trop=zeros(size(Flux_model_lut_depth));
for i=idx_open(idx_trop)
    in = inpolygon(lon_flux_lut_180,lat_flux_lut,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_lutz_trop(in==1)=1;
end
idx_lutz_trop=logical(idx_lutz_trop);


% figure
% scatter(lon_flux_lut_180,lat_flux_lut,8,'b','filled');
% hold on
% scatter(lon_flux_lut_180(idx_lutz_trop),lat_flux_lut(idx_lutz_trop),8,'r','filled');
% scatter(lon_flux_lut_180(idx_lutz_oligo),lat_flux_lut(idx_lutz_oligo),8,'r','filled');
% scatter(lon_flux_lut_180(idx_lutz_temp),lat_flux_lut(idx_lutz_temp),8,'r','filled');
% scatter(lon_flux_lut_180(idx_lutz_polar),lat_flux_lut(idx_lutz_polar),8,'r','filled');
% geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
% title('Carbon export data locations compiled in Lutz et al 2007')
% set(gcf,'color','w');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CARBON EXPORT - DATA FROM HENSON 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat=importdata('data_carbon_export/gbc20895-sup-0002-2018gb006158-data set s1.xlsx');

lat_data_flux=dat.data(:,1);
lon_data_flux=dat.data(:,2);
day_data_flux=dat.data(:,3);
flux_data=dat.data(:,4);
PP_data_flux=dat.data(:,5);
eratio_data=dat.data(:,6);
depth_data=dat.data(:,7);
macrozoo_data_flux=dat.data(:,15);
mesozoo_data_flux=dat.data(:,13);

z_interp_s=0:240;
zface_interp=zface;

lon2_flux=mod(lon_data_flux(:),360);

zface=cumsum(deltaz);

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
monticks = [(1+cumsum(mon)) 365];

NPP_t_100=squeeze(nansum(NPP_t(:,:,1:2,:),3));

flux_model=NaN(length(lon_data_flux),1);
Flux_model_hens=NaN(length(lon_data_flux),1);
eff_model=NaN(length(lon_data_flux),1);
PP_model_flux=NaN(length(lon_data_flux),1);
for i=1:length(lon_data_flux)
    
    if isnan(lat_data_flux(i)) || isnan(lat_data_flux(i))
        xid=find(x>lon2_flux(i) & x<lon2_flux(i)+2.9);
        yid=find(y>lat_data_flux(i) & y<lat_data_flux(i)+2.9);
        flux_model(i)=NaN;
        eff_model(i)=NaN;
        PP_model_flux(i)=NaN;
    else
        xid=find(x>lon2_flux(i) & x<lon2_flux(i)+2.9);
        if lon2_flux(i)>max(x)
            xid=128;
        end
        yid=find(y>lat_data_flux(i) & y<lat_data_flux(i)+2.9);
        if depth_data(i)<130
            flux_model(i)=nanmean(Export_t(xid(1),yid(1),2,day_data_flux(i)),4);           
        elseif depth_data(i)>=130
            flux_model(i)=nanmean(Export_t(xid(1),yid(1),3,day_data_flux(i)),4);
        end
        F_wc=[squeeze(nanmean(Export_t(xid(1),yid(1),:,day_data_flux(i)),4))];
        if ~isnan(F_wc(2))
            F_interp = interp1(zface_interp,F_wc,z_interp_s);  
            idx_z=find(z_interp_s>depth_data(i));
            Flux_model_hens(i)=F_interp(idx_z(1));
%             Flux_model_hens_depth_interp(i,:)=F_interp;
        end
        
        eff_model(i)=flux_model(i)./nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
        PP_model_flux(i)=nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
    end

    
end


idx_data=find(~isnan(eff_model) & ~isinf(eff_model) & ~isnan(eratio_data) & PP_model_flux>=0 & flux_data>0);

flux_model2=flux_model(idx_data);
flux_model_hens2=Flux_model_hens(idx_data);
eff_model2=eff_model(idx_data);
PP_model_flux2=PP_model_flux(idx_data);
flux_data2=flux_data(idx_data);
eratio2=eratio_data(idx_data);
lat_data2=lat_data_flux(idx_data);
lon_data2=lon_data_flux(idx_data);
PP_data_flux2=PP_data_flux(idx_data);



%% get Ocean basins
%1- polar north
% 2:21 - Atlantic ocean
% 22:30 - Indian ocean
% 31:50 - Pacific ocean
% 50:54 - Southern ocean

s = shaperead('../3Doutputs/Longhurst/longhurst_v4_2010/Longhurst_world_v4_2010.shp');

idx_polar=[32,31,3,2,1]; 
idx_temp=[30,29,20,18,17,16,13,5,4]; 
idx_trop=[26,25,14,9,8,27]; 
idx_oligo=[24,23,21,15,10,7,6];

%remove coastal regions
idx_open=zeros(1,54);
for i=1:54
   alk=split(s(i).ProvDescr);
   if ~contains(alk{1},'Coastal')
      idx_open(i)=i;
   end   
end
idx_open=idx_open(idx_open~=0);

idx_hens_oligo=zeros(size(Flux_model_hens));
for i=idx_open(idx_oligo)
    in = inpolygon(lon_data2,lat_data2,s(i).X,s(i).Y);
    in=in';
    idx_hens_oligo(in==1)=1;
end
idx_hens_oligo=logical(idx_hens_oligo);

idx_hens_polar=zeros(size(Flux_model_hens));
for i=idx_open(idx_polar)
    in = inpolygon(lon_data2,lat_data2,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_hens_polar(in==1)=1;
end
idx_hens_polar=logical(idx_hens_polar);


idx_hens_temp=zeros(size(Flux_model_hens));
for i=idx_open(idx_temp)
    in = inpolygon(lon_data2,lat_data2,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_hens_temp(in==1)=1;
end
idx_hens_temp=logical(idx_hens_temp);

idx_hens_trop=zeros(size(Flux_model_hens));
for i=idx_open(idx_trop)
    in = inpolygon(lon_data2,lat_data2,s(i).X,s(i).Y);
    in=in';
    % returns in indicating if the query points specified by xq and yq are inside or on the edge of the polygon area defined by xv and yv.
    idx_hens_trop(in==1)=1;
end
idx_hens_trop=logical(idx_hens_trop);

idx_hens_used=idx_hens_trop | idx_hens_oligo | idx_hens_temp | idx_hens_polar;


%% binned carbon export from henson

flux_data_tmm=NaN(length(x),length(y),365);
flux_data_tmm_std=NaN(length(x),length(y),365);
flux_depth_tmm=NaN(length(x),length(y),365);
flux_lon_tmm=NaN(length(x),length(y),365);
flux_lat_tmm=NaN(length(x),length(y),365);
for i=1:length(x)-1
    for j=1:length(y)-1
        int_x=[x(i), x(i+1)];
        int_y=[y(j), y(j+1)];
        
        for k=1:365
            idx_lonlat=find(lon2_flux>=int_x(1) & lon2_flux<int_x(2) & lat_data_flux>=int_y(1) & lat_data_flux<int_y(2) & day_data_flux==k);
                if ~isempty(idx_lonlat)
                    flux_data_tmm(i,j,k)=nanmean(flux_data(idx_lonlat));
                    flux_data_tmm_std(i,j,k)=nanstd(flux_data(idx_lonlat));
                    flux_depth_tmm(i,j,k)=nanmean(depth_data(idx_lonlat));
                    flux_lon_tmm(i,j,k)=x(i);
                    flux_lat_tmm(i,j,k)=y(j);
                end
        end
        
    end
end

flux_model=squeeze(Export_t(:,:,2,:));
flux_model3=squeeze(Export_t(:,:,3,:));
flux_model(flux_depth_tmm>130)=flux_model3(flux_depth_tmm>130);

load('tempday.mat')
sst_t=matrixToGrid(tempday(:,1:2:end), [], 'boxes.mat', 'grid.mat');
sst_3d=squeeze(mean(sst_t(:,:,1,:),3));

corr([flux_data_tmm(:),flux_model(:)],'Rows','complete')

corr([flux_data_tmm(sst_3d>5),flux_model(sst_3d>5)],'Rows','complete')


flux_data_4=flux_data_tmm;%(flux_data_tmm_std~=0);
flux_model_4=flux_model;%(flux_data_tmm_std~=0);
flux_data_4(flux_data_4<0)=NaN;
flux_lon_4=flux_lon_tmm;%(flux_data_tmm_std~=0);
flux_lat_4=flux_lat_tmm;%(flux_data_tmm_std~=0);

flux_data_5=flux_data_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
flux_model_5=flux_model_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
flux_lon_5=flux_lon_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
flux_lat_5=flux_lat_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));

figure
subplot(2,1,1)
scatter(lon_data2,lat_data2,8,'r','filled');
hold on
scatter(lon_data2(idx_hens_used),lat_data2(idx_hens_used),8,'b','filled');
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
title('Carbon export data locations compiled in Henson et al 2019')
set(gcf,'color','w');

subplot(2,1,2)
scatter(lon_flux_lut_180,lat_flux_lut,8,'r','filled');
hold on
scatter(lon_flux_lut_180(idx_lutz_trop),lat_flux_lut(idx_lutz_trop),8,'b','filled');
scatter(lon_flux_lut_180(idx_lutz_oligo),lat_flux_lut(idx_lutz_oligo),8,'b','filled');
scatter(lon_flux_lut_180(idx_lutz_temp),lat_flux_lut(idx_lutz_temp),8,'b','filled');
scatter(lon_flux_lut_180(idx_lutz_polar),lat_flux_lut(idx_lutz_polar),8,'b','filled');
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
title('Carbon export data locations compiled in Lutz et al 2007')
set(gcf,'color','w');


% print -depsc figures_paper/fig3_SI_export_data

%% Plots data comparison NPP and export

st=1;
x0=0;
y0=0;
width=18;
height=10;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');


tiledlayout(2,4,'Padding','Compact','TileSpacing','Compact')

nexttile
ssiz=3;
cmap=brewermap(6,'Dark2');
% errorbar(NPP_tmm(idx_SO),NPP_int(idx_SO),NPP_tmm_std(idx_SO),'horizontal','o'...
%     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% hold on
% errorbar(NPP_tmm(idx_HOT),NPP_int(idx_HOT),NPP_tmm_std(idx_HOT),'horizontal','o'...
%     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% errorbar(NPP_tmm(idx_BATS),NPP_int(idx_BATS),NPP_tmm_std(idx_BATS),'horizontal','o'...
%     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% errorbar(NPP_tmm(idx_NA),NPP_int(idx_NA),NPP_tmm_std(idx_NA),'horizontal','o'...
%     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% errorbar(NPP_tmm(idx_arab),NPP_int(idx_arab),NPP_tmm_std(idx_arab),'horizontal','o'...
%     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
hold on
h1=scatter(NPP_tmm(idx_SO),NPP_int(idx_SO),ssiz,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5);
h2=scatter(NPP_tmm(idx_HOT),NPP_int(idx_HOT),ssiz,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5);
h3=scatter(NPP_tmm(idx_BATS),NPP_int(idx_BATS),ssiz,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5);
h4=scatter(NPP_tmm(idx_NA),NPP_int(idx_NA),ssiz,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5);
h5=scatter(NPP_tmm(idx_arab),NPP_int(idx_arab),ssiz,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:),'markerfacealpha',0.5);
plot(linspace(0,2000),linspace(0,2000),'k--')
xlabel('Data [mgC m^{-2} day^{-1}]')
ylabel('Model [mgC m^{-2} day^{-1}]')
title('{\bfa.} NPP','fontweight','normal')
corval=corr(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)),'type','Spearman','Rows','complete');
% text(0.05,0.95,strcat(['r=',num2str(corval,'%.2f')]),'Units','normalized','fontsize',9)
rmse= fitlm(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)));
Err= (1/length(~isnan(NPP_tmm)))*(nansum((NPP_int(:))-(NPP_tmm(:))));
ax = gca;
% ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e2 1.3e3])
xlim([5e1 5e3])
yticks([1e2 1e3])
xticks([1e2 1e3])
leg1=legend([h1,h2,h3,h4,h5],{'SO','HOT','BATS','NA','Arab'});
% leg1=legend([h1,h2,h3,h4,h5],{'SO','HOT','BATS','NA','Arab'});
% leg1.Title.String=('NPP');
title(leg1,'NPP')
legend boxoff
leg1.Layout.Tile = 4;
% % leg1.ItemTokenSize = [10,5]; 
% icons = findobj(icons,'Type','patch');
% icons = findobj(icons,'Marker','none','-xor');
% set(icons,'MarkerSize',5);
% % title(leg1,'NPP')
% htitle = get(leg1,'Title');
% set(htitle,'String','NPP')
set(gcf,'color','w');

idx=isfinite(NPP_tmm(idx_HOT)) & isfinite(NPP_int(idx_HOT));
corrcoef(log10(NPP_tmm(idx_HOT(idx))),log10(NPP_int(idx_HOT(idx))))
fitlm(log10(NPP_tmm(idx_HOT(idx))),log10(NPP_int(idx_HOT(idx))))
[~,~,~, logMAE]=func_calc_bias(NPP_int(idx_HOT(idx)),NPP_tmm(idx_HOT(idx)))

idx_loc=idx_NA;
idx=isfinite(NPP_tmm(idx_loc)) & isfinite(NPP_int(idx_loc));
corrcoef(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
fitlm(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
[~,~,~, logMAE]=func_calc_bias(NPP_int(idx_loc(idx)),NPP_tmm(idx_loc(idx)))

idx_loc=idx_BATS;
idx=isfinite(NPP_tmm(idx_loc)) & isfinite(NPP_int(idx_loc));
corrcoef(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
fitlm(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
[~,~,~, logMAE]=func_calc_bias(NPP_int(idx_loc(idx)),NPP_tmm(idx_loc(idx)))

idx_loc=idx_arab;
idx=isfinite(NPP_tmm(idx_loc)) & isfinite(NPP_int(idx_loc));
corrcoef(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
fitlm(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
[~,~,~, logMAE]=func_calc_bias(NPP_int(idx_loc(idx)),NPP_tmm(idx_loc(idx)))

idx_loc=idx_SO;
idx=isfinite(NPP_tmm(idx_loc)) & isfinite(NPP_int(idx_loc));
corrcoef(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
fitlm(log10(NPP_tmm(idx_loc(idx))),log10(NPP_int(idx_loc(idx))))
[~,~,~, logMAE]=func_calc_bias(NPP_int(idx_loc(idx)),NPP_tmm(idx_loc(idx)))




idx=isfinite(flux_data2(idx_hens_trop)) & isfinite(flux_model_hens2(idx_hens_trop));
corrcoef(log10(flux_data2(idx_hens_trop)),log10(flux_model_hens2(idx_hens_trop)))
fitlm(log10(flux_data2(idx_hens_trop)),log10(flux_model_hens2(idx_hens_trop)))
[~,~,~, logMAE]=func_calc_bias(flux_model_hens2(idx_hens_trop),flux_data2(idx_hens_trop))

Modl=flux_model_hens2(idx_hens_oligo);
Obs=flux_data2(idx_hens_oligo);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl),log10(Obs))
fitlm(log10(Modl),log10(Obs))
[~,~,~, logMAE]=func_calc_bias(Obs,Modl)

Modl=flux_model_hens2(idx_hens_temp);
Obs=flux_data2(idx_hens_temp);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

Modl=flux_model_hens2(idx_hens_polar);
Obs=flux_data2(idx_hens_polar);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

cmap=brewermap(4,'Spectral');

nexttile
scatter(flux_data2(idx_hens_trop),flux_model_hens2(idx_hens_trop),ssiz,'filled','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
scatter(flux_data2(idx_hens_oligo),flux_model_hens2(idx_hens_oligo),ssiz,'filled','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
scatter(flux_data2(idx_hens_temp),flux_model_hens2(idx_hens_temp),ssiz,'filled','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
scatter(flux_data2(idx_hens_polar),flux_model_hens2(idx_hens_polar),ssiz,'filled','markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5)

% scatter(flux_data2(idx_hens_used),flux_model_hens2(idx_hens_used),ssiz,lat_data2(idx_hens_used))
% hold on
plot(linspace(1e0,1e3),linspace(1e0,1e3),'k--')
xlabel('Data [mgC m^{-2} day^{-1}]')
set(gca,'yscale','log')
set(gca,'xscale','log')
% colorbar
title('Carbon export')
ylim([5 1e3])
xlim([5 1e3])
xticks([1e1 1e2 1e3])
yticks([1e1 1e2 1e3])
title('{\bfb.} Surface export (^{243}Th)','fontweight','normal')


Modl=Flux_model_lut(idx_lutz_polar);
Obs=flux_data_lut(idx_lutz_polar);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

Modl=Flux_model_lut(idx_lutz_oligo);
Obs=flux_data_lut(idx_lutz_oligo);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

Modl=Flux_model_lut(idx_lutz_trop);
Obs=flux_data_lut(idx_lutz_trop);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

Modl=Flux_model_lut(idx_lutz_temp);
Obs=flux_data_lut(idx_lutz_temp);
idx=isfinite(Modl) & isfinite(Obs);
corrcoef(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

Modl=Flux_model_lut;
Obs=flux_data_lut;
idx=isfinite(Modl) & isfinite(Obs);
corr(log10(Modl(idx)),log10(Obs(idx)))
fitlm(log10(Modl(idx)),log10(Obs(idx)))
[~,~,~, logMAE]=func_calc_bias(Obs(idx),Modl(idx))

nexttile
scatter(flux_data_lut(idx_lutz_trop),Flux_model_lut(idx_lutz_trop),ssiz,'filled','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
scatter(flux_data_lut(idx_lutz_oligo),Flux_model_lut(idx_lutz_oligo),ssiz,'filled','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
scatter(flux_data_lut(idx_lutz_temp),Flux_model_lut(idx_lutz_temp),ssiz,'filled','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
scatter(flux_data_lut(idx_lutz_polar),Flux_model_lut(idx_lutz_polar),ssiz,'filled','markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5)
plot(logspace(-1,2),logspace(-1,2),'k--')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Data [mgC m^{-2} day^{-1}]')
title('{\bfc.} Deep export (traps)','fontweight','normal')
leg2=legend('Tropical','Oligotrophic','Temperate','Polar');
title(leg2,'Export')
legend boxoff
% icons = findobj(icons,'Type','patch');
% icons = findobj(icons,'Marker','none','-xor');
% set(icons,'MarkerSize',5);
% leg2.Layout.Tile = 4;
% ylim([5e-1 5e1])
xticks([1e-1 1e0 1e1 1e2])
yticks([1e-1 1e0 1e1 1e2])
[~,~,~, logMAE]=func_calc_bias(Flux_model_lut(idx_lutz_polar),flux_data_lut(idx_lutz_polar))


nexttile(5)
% scatter(flux_data_lut(idx_lutz_trop),-depth_flux_lut(idx_lutz_trop),5,'k','filled')
hold on
% plot(nanmean(Flux_model_lut_depth_interp(idx_lutz_trop,:),1),-z_interp','r')
curve1=nanmean(Flux_model_lut_depth_interp(idx_lutz_trop,:),1)-nanstd(Flux_model_lut_depth_interp(idx_lutz_trop,:),1);
curve2=nanmean(Flux_model_lut_depth_interp(idx_lutz_trop,:),1)+nanstd(Flux_model_lut_depth_interp(idx_lutz_trop,:),1);
idx=isfinite(curve1) & isfinite(curve2) & curve1>0;
cpl=patch([curve1(idx), fliplr(curve2(idx))],[-z_interp(idx), fliplr(-z_interp(idx))],'k','EdgeColor','none','facealpha',0.1);
scatter(flux_data_lut(idx_lutz_trop),-depth_flux_lut(idx_lutz_trop),5,'k','filled','markerfacecolor',cmap(1,:),'markeredgecolor',[92, 92, 92]./255)
hold on
plot(nanmean(Flux_model_lut_depth_interp(idx_lutz_trop,:),1),-z_interp','k','linewidth',1)
% plot(nanmean(Export_y_trop2,1),-zface','k--','linewidth',1)
set(gca,'xscale','log')
xlim([10^-0.5 10^2.5])
% xlabel('Export [mgC m^{-2} d^{-1}]')
ylabel('Depth [Km]')
title('Tropical')
title('{\bfd.} Tropical','fontweight','normal')
ylim([-6000 0])
yticks([-6e3 -4e3 -2e3 0])
yticklabels({'-6','-4','-2','0'})
xticks([1e0 1e1 1e2])
xlabel('Export [mgC m^{-2} d^{-1}]')



nexttile(6)
% scatter(flux_data_lut(idx_lutz_oligo),-depth_flux_lut(idx_lutz_oligo),5,'k','filled')
hold on
% plot(Flux_model_lut_depth(idx_lutz_oligo,:),-zface','r')
curve1=nanmean(Flux_model_lut_depth_interp(idx_lutz_oligo,:),1)-nanstd(Flux_model_lut_depth_interp(idx_lutz_oligo,:),1);
curve2=nanmean(Flux_model_lut_depth_interp(idx_lutz_oligo,:),1)+nanstd(Flux_model_lut_depth_interp(idx_lutz_oligo,:),1);
idx=isfinite(curve1) & isfinite(curve2) & curve1>0;
cpl=patch([curve1(idx), fliplr(curve2(idx))],[-z_interp(idx), fliplr(-z_interp(idx))],'k','EdgeColor','none','facealpha',0.1);
hold on
scatter(flux_data_lut(idx_lutz_oligo),-depth_flux_lut(idx_lutz_oligo),5,'k','filled','markerfacecolor',cmap(2,:),'markeredgecolor',[92, 92, 92]./255)
plot(nanmean(Flux_model_lut_depth_interp(idx_lutz_oligo,:),1),-z_interp','k','linewidth',1)
% plot(nanmean(Export_y_oligo2,1),-zface','k--','linewidth',1)
set(gca,'xscale','log')
xlim([10^-0.5 10^2.5])
xlabel('Export [mgC m^{-2} d^{-1}]')
% title('Oligotrophic')
title('{\bfe.} Oligotrohic','fontweight','normal')
ylim([-6000 0])
yticks([-6e3 -4e3 -2e3 0])
yticklabels({'-6','-4','-2','0'})
xticks([1e0 1e1 1e2])



nexttile(7)
% scatter(flux_data_lut(idx_lutz_temp),-depth_flux_lut(idx_lutz_temp),5,'k','filled')
hold on
% plot(Flux_model_lut_depth(idx_lutz_temp,:),-zface','r')
curve1=nanmean(Flux_model_lut_depth_interp(idx_lutz_temp,:),1)-nanstd(Flux_model_lut_depth_interp(idx_lutz_temp,:),1);
curve2=nanmean(Flux_model_lut_depth_interp(idx_lutz_temp,:),1)+nanstd(Flux_model_lut_depth_interp(idx_lutz_temp,:),1);
idx=isfinite(curve1) & isfinite(curve2) & curve1>0;
cpl=patch([curve1(idx), fliplr(curve2(idx))],[-z_interp(idx), fliplr(-z_interp(idx))],'k','EdgeColor','none','facealpha',0.1);
hold on
scatter(flux_data_lut(idx_lutz_temp),-depth_flux_lut(idx_lutz_temp),5,'k','filled','markerfacecolor',cmap(3,:),'markeredgecolor',[92, 92, 92]./255)
plot(nanmean(Flux_model_lut_depth_interp(idx_lutz_temp,:),1),-z_interp','k','linewidth',1)
% plot(nanmean(Export_y_temp2,1),-zface','k--','linewidth',1)
set(gca,'xscale','log')
xlim([10^-0.5 10^2.5])
% xlabel('Export [mgC m^{-2} d^{-1}]')
title('{\bff.} Temperate','fontweight','normal')
% yticks([])
% yticklabels([])
ylim([-6000 0])
yticks([-6e3 -4e3 -2e3 0])
yticklabels({'-6','-4','-2','0'})
xticks([1e0 1e1 1e2])
xlabel('Export [mgC m^{-2} d^{-1}]')




nexttile(8)
hold on
% plot(Flux_model_lut_depth(idx_lutz_polar,:),-zface','r')
curve1=nanmean(Flux_model_lut_depth_interp(idx_lutz_polar,:),1)-nanstd(Flux_model_lut_depth_interp(idx_lutz_polar,:),1);
curve2=nanmean(Flux_model_lut_depth_interp(idx_lutz_polar,:),1)+nanstd(Flux_model_lut_depth_interp(idx_lutz_polar,:),1);
idx=isfinite(curve1) & isfinite(curve2) & curve1>0;
cpl=patch([curve1(idx), fliplr(curve2(idx))],[-z_interp(idx), fliplr(-z_interp(idx))],'k','EdgeColor','none','facealpha',0.1);
% camroll(180)
% plot([curve1;curve2], [-z_interp;-z_interp],'color',[217, 219, 218]./244);
hold on
scatter(flux_data_lut(idx_lutz_polar),-depth_flux_lut(idx_lutz_polar),5,'k','filled','markerfacecolor',cmap(4,:),'markeredgecolor',[92, 92, 92]./255)
plot(nanmean(Flux_model_lut_depth_interp(idx_lutz_polar,:),1),-z_interp','k','linewidth',1)
% plot(nanmean(Export_y_polar2,1),-zface','k--','linewidth',1)
set(gca,'xscale','log')
xlim([10^-0.5 10^2.5])
% xlabel('Export [mgC m^{-2} d^{-1}]')
title('{\bfg.} Polar','fontweight','normal')
ylim([-6000 0])
yticks([-6e3 -4e3 -2e3 0])
yticklabels({'-6','-4','-2','0'})
xticks([1e0 1e1 1e2])
xlabel('Export [mgC m^{-2} d^{-1}]')

set(gcf,'color','w')

% print -depsc figures_paper/fig3_NPP_data
print(gcf, 'fig3_nppdat.pdf', '-dpdf', '-fillpage')

%% NPP all data supplementary information

st=1;
x0=0;
y0=0;
width=8;
height=7;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

ssiz=3;
cmap=brewermap(6,'Dark2');
hold on
h1=scatter(NPP_tmm(idx_SO),NPP_int(idx_SO),ssiz,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5);
h2=scatter(NPP_tmm(idx_HOT),NPP_int(idx_HOT),ssiz,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5);
h3=scatter(NPP_tmm(idx_BATS),NPP_int(idx_BATS),ssiz,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5);
h4=scatter(NPP_tmm(idx_NA),NPP_int(idx_NA),ssiz,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5);
h5=scatter(NPP_tmm(idx_arab),NPP_int(idx_arab),ssiz,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:),'markerfacealpha',0.5);
plot(linspace(1e0,1e4),linspace(1e0,1e4),'--','Color','k')
xlabel('Data [mgC m^{-2} day^{-1}]')
ylabel('Model [mgC m^{-2} day^{-1}]')
title('NPP')
corval=corr(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)),'type','Spearman','Rows','complete');
% text(0.05,0.95,strcat(['r=',num2str(corval,'%.2f')]),'Units','normalized','fontsize',9)
rmse= fitlm(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)));
Err= (1/length(~isnan(NPP_tmm)))*(nansum((NPP_int(:))-(NPP_tmm(:))));
ax = gca;
% ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e0 1e4])
xlim([1e0 1e4])
% yticks([2e2 1e3])
leg1=legend([h1,h2,h3,h4,h5],{'SO','HOT','BATS','NA','Arab'},'Location','northwest');
% leg1=legend([h1,h2,h3,h4,h5],{'SO','HOT','BATS','NA','Arab'});
% leg1.Title.String=('NPP');
% title(leg1,'NPP')
legend boxoff

set(gcf,'color','w');
% print -depsc figures_paper/fig3_NPP_data_SO

%%
% % %% binned carbon export
% % % several data points fall within the same TM bin and in the same day. We
% % % combinethose and do the mean
% % flux_data_tmm=NaN(length(x),length(y),365);
% % flux_data_tmm_std=NaN(length(x),length(y),365);
% % flux_depth_tmm=NaN(length(x),length(y),365);
% % flux_lon_tmm=NaN(length(x),length(y),365);
% % flux_lat_tmm=NaN(length(x),length(y),365);
% % for i=1:length(x)-1
% %     for j=1:length(y)-1
% %         int_x=[x(i), x(i+1)];
% %         int_y=[y(j), y(j+1)];
% %         
% %         for k=1:365
% %             idx_lonlat=find(lon2_flux>=int_x(1) & lon2_flux<int_x(2) & lat_data_flux>=int_y(1) & lat_data_flux<int_y(2) & day_data_flux==k);
% %                 if ~isempty(idx_lonlat)
% %                     flux_data_tmm(i,j,k)=nanmean(flux_data(idx_lonlat));
% %                     flux_data_tmm_std(i,j,k)=nanstd(flux_data(idx_lonlat));
% %                     flux_depth_tmm(i,j,k)=nanmean(depth_data(idx_lonlat));
% %                     flux_lon_tmm(i,j,k)=x(i);
% %                     flux_lat_tmm(i,j,k)=y(j);
% %                 end
% %         end
% %         
% %     end
% % end
% % 
% % flux_model=squeeze(Export_t(:,:,2,:));
% % flux_model3=squeeze(Export_t(:,:,3,:));
% % flux_model(flux_depth_tmm>130)=flux_model3(flux_depth_tmm>130);
% % 
% % load('tempday.mat')
% % sst_t=matrixToGrid(tempday(:,1:2:end), [], 'boxes.mat', 'grid.mat');
% % sst_3d=squeeze(mean(sst_t(:,:,1,:),3));
% % 
% % corr([flux_data_tmm(:),flux_model(:)],'Rows','complete')
% % 
% % corr([flux_data_tmm(sst_3d>5),flux_model(sst_3d>5)],'Rows','complete')
% % 
% % 
% % flux_data_4=flux_data_tmm;%(flux_data_tmm_std~=0);
% % flux_model_4=flux_model;%(flux_data_tmm_std~=0);
% % flux_data_4(flux_data_4<0)=NaN;
% % flux_lon_4=flux_lon_tmm;%(flux_data_tmm_std~=0);
% % flux_lat_4=flux_lat_tmm;%(flux_data_tmm_std~=0);
% % 
% % flux_data_5=flux_data_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
% % flux_model_5=flux_model_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
% % flux_lon_5=flux_lon_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
% % flux_lat_5=flux_lat_4(~isnan(flux_data_4) & ~isnan(flux_model_4) & ~isinf(log10(flux_data_4)));
% % 
% % figure
% % scatter(lon_data2,lat_data2,8,'r','filled');
% % hold on
% % geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
% % title('Carbon export data location')
% % set(gcf,'color','w');
% % 
% % 
% % %% with binned data
% % 
% % st=1;
% % x0=0;
% % y0=0;
% % width=18;
% % height=7;
% % fig=figure('Units','centimeters',...
% % 'Position',[x0 y0 width height],...
% % 'PaperPositionMode','auto');
% % 
% % greycolor=[0.5 0.5 0.5]
% % 
% % subplot(1,2,1)
% % ssiz=10;
% % cmap=brewermap(6,'Dark2');
% % errorbar(NPP_tmm(idx_SO),NPP_int(idx_SO),NPP_tmm_std(idx_SO),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % hold on
% % errorbar(NPP_tmm(idx_HOT),NPP_int(idx_HOT),NPP_tmm_std(idx_HOT),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_BATS),NPP_int(idx_BATS),NPP_tmm_std(idx_BATS),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_NA),NPP_int(idx_NA),NPP_tmm_std(idx_NA),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_arab),NPP_int(idx_arab),NPP_tmm_std(idx_arab),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % h1=scatter(NPP_tmm(idx_SO),NPP_int(idx_SO),ssiz,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5);
% % h2=scatter(NPP_tmm(idx_HOT),NPP_int(idx_HOT),ssiz,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5);
% % h3=scatter(NPP_tmm(idx_BATS),NPP_int(idx_BATS),ssiz,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5);
% % h4=scatter(NPP_tmm(idx_NA),NPP_int(idx_NA),ssiz,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5);
% % h5=scatter(NPP_tmm(idx_arab),NPP_int(idx_arab),ssiz,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:),'markerfacealpha',0.5);
% % plot(linspace(0,2000),linspace(0,2000),'--','Color',[0.5 0.5 0.5])
% % legend([h1,h2,h3,h4,h5],'SO','HOT','BATS','NA','Arab')
% % set(gcf,'color','w');
% % xlabel('Data [mgC m^{-2} day^{-1}]','fontsize',10)
% % ylabel('Model [mgC m^{-2} day^{-1}]','fontsize',10)
% % title('{\bfa.} NPP by regions','fontweight','normal')
% % corval=corr(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)),'type','Spearman','Rows','complete');
% % % text(0.05,0.95,strcat(['r=',num2str(corval,'%.2f')]),'Units','normalized','fontsize',9)
% % rmse= fitlm(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)));
% % Err= (1/length(~isnan(NPP_tmm)))*(nansum((NPP_int(:))-(NPP_tmm(:))));
% % ax = gca;
% % % ax.TitleHorizontalAlignment = 'left';
% % set(gca,'box','off') 
% % 
% % 
% % subplot(1,2,2)
% % errorbar(flux_data_tmm(flux_data_tmm_std~=0),flux_model(flux_data_tmm_std~=0),flux_data_tmm_std(flux_data_tmm_std~=0),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % hold on
% % scatter(flux_data_tmm(:),flux_model(:),10,sst_3d(:),'filled','markerfacealpha',0.5)
% % plot(logspace(0,3),logspace(0,3),'--','color',[0.5 0.5 0.5])
% % colormap(viridis(100))
% % ch=colorbar;
% % ylabel(ch,'Latitude')
% % set(gca,'yscale','log')
% % set(gca,'xscale','log')
% % ylim([10^0.5 1e3])
% % xlim([10^0.5 1e3])
% % set(gcf,'color','w');
% % xlabel('Data [mgC m^{-2} day^{-2}]')
% % ylabel('Model [mgC m^{-2} day^{-2}]')
% % title('{\bfb.} Carbon export','fontweight','normal')
% % ax = gca;
% % set(gca,'box','off') 
% % 
% % 
% % % print -depsc figures_paper/fig3_NPP_data
% % 
% % %%
% % 
% % %% with binned data
% % 
% % st=1;
% % x0=0;
% % y0=0;
% % width=18;
% % height=7;
% % fig=figure('Units','centimeters',...
% % 'Position',[x0 y0 width height],...
% % 'PaperPositionMode','auto');
% % 
% % greycolor=[0.5 0.5 0.5]
% % 
% % subplot(1,2,1)
% % ssiz=10;
% % cmap=brewermap(6,'Dark2');
% % errorbar(NPP_tmm(idx_SO),NPP_int(idx_SO),NPP_tmm_std(idx_SO),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % hold on
% % errorbar(NPP_tmm(idx_HOT),NPP_int(idx_HOT),NPP_tmm_std(idx_HOT),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_BATS),NPP_int(idx_BATS),NPP_tmm_std(idx_BATS),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_NA),NPP_int(idx_NA),NPP_tmm_std(idx_NA),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % errorbar(NPP_tmm(idx_arab),NPP_int(idx_arab),NPP_tmm_std(idx_arab),'horizontal','o'...
% %     ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
% % h1=scatter(NPP_tmm(idx_SO),NPP_int(idx_SO),ssiz,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5);
% % h2=scatter(NPP_tmm(idx_HOT),NPP_int(idx_HOT),ssiz,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5);
% % h3=scatter(NPP_tmm(idx_BATS),NPP_int(idx_BATS),ssiz,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5);
% % h4=scatter(NPP_tmm(idx_NA),NPP_int(idx_NA),ssiz,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5);
% % h5=scatter(NPP_tmm(idx_arab),NPP_int(idx_arab),ssiz,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:),'markerfacealpha',0.5);
% % plot(linspace(0,2000),linspace(0,2000),'--','Color',[0.5 0.5 0.5])
% % legend([h1,h2,h3,h4,h5],'SO','HOT','BATS','NA','Arab')
% % set(gcf,'color','w');
% % xlabel('Data [mgC m^{-2} day^{-1}]','fontsize',10)
% % ylabel('Model [mgC m^{-2} day^{-1}]','fontsize',10)
% % title('{\bfa.} NPP by regions','fontweight','normal')
% % corval=corr(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)),'type','Spearman','Rows','complete');
% % % text(0.05,0.95,strcat(['r=',num2str(corval,'%.2f')]),'Units','normalized','fontsize',9)
% % rmse= fitlm(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)));
% % Err= (1/length(~isnan(NPP_tmm)))*(nansum((NPP_int(:))-(NPP_tmm(:))));
% % ax = gca;
% % % ax.TitleHorizontalAlignment = 'left';
% % set(gca,'box','off') 
% % 
% % 
% % iidx_SO= isfinite(NPP_int) & isfinite(NPP_tmm) & (lat_3D<-50);
% % iidx_HOT= isfinite(NPP_int) & isfinite(NPP_tmm) & (lon_3D<-150 & lat_3D>0);
% % iidx_BATS= isfinite(NPP_int) & isfinite(NPP_tmm) & (lon_3D>-67 & lon_3D<-65 & lat_3D>29 & lat_3D<30);
% % iidx_NA= isfinite(NPP_int) & isfinite(NPP_tmm) & (lon_3D>-21 & lon_3D<0);
% % iidx_arab= isfinite(NPP_int) & isfinite(NPP_tmm) & (lon_3D>50 & lon_3D<67);
% % 
% % aldat=[geomean(NPP_tmm(iidx_arab)), geomean(NPP_tmm(iidx_HOT)), geomean(NPP_tmm(iidx_BATS)), geomean(NPP_tmm(iidx_NA)), geomean(NPP_tmm(iidx_SO))];
% % almod=[geomean(NPP_int(iidx_arab)), geomean(NPP_int(iidx_HOT)), geomean(NPP_int(iidx_BATS)), geomean(NPP_int(iidx_NA)), geomean(NPP_int(iidx_SO))];
% % 
% % figure
% % scatter(almod,aldat)
% % hold on
% % plot(linspace(0,1200),linspace(0,1200))
% % 
% % idx= isfinite(NPP_int) & isfinite(NPP_tmm) & (iidx_arab | iidx_NA | iidx_BATS | iidx_HOT);% | idx_SO;
% % 
% % figure
% % scatter(NPP_int(idx),NPP_tmm(idx))
% % 
% % 
% % 
% % figure
% % scatter(log10(NPP_tmm(:)),log10(NPP_int(:)))

%% rates compared to data plots

% addpath('C:/Users/Camila/Desktop/Backup/Projects/Cluster runs/cmap_ocean/')

x0=0;
y0=0;
width=17;
height=11;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
% tiledlayout(2,6,'Padding','Compact');

nexttile
% subplot(2,3,1)
h=function_plot_map_contours(nansum(nansum(NPP_t,3),4)./1000,0:50:500);
ylabel(h,'[gC m^{-2} year^{-1}]')
% caxis([0 9])
title('{\bfa.} NPP','fontweight','normal','fontsize',10)
caxis([0 500])

nexttile
% subplot(2,3,2)
h=function_plot_map_contours(nansum(Export_t(:,:,2,:),4)./1000,0:5:50);
ylabel(h,'[gC m^{-2} year^{-1}]')
title('{\bfb.} Export 120m','fontweight','normal','fontsize',10)
caxis([0 50])

nexttile
% subplot(2,3,3)
h=function_plot_map_contours(nansum(Export_t(:,:,7,:),4)./1000,0:1:10);
ylabel(h,'[gC m^{-2} year^{-1}]')
title('{\bfc.} Export 1080m','fontweight','normal','fontsize',10)
caxis([0 10])

nexttile(5)
% subplot(2,3,5)
 pe=nansum(Export_t(:,:,2,:),4)./nansum(nansum(NPP_t,3),4);
  pe(isinf(pe))=NaN;
h=function_plot_map_contours(pe,0:0.04:0.4);
caxis([0 0.4])
ylabel(h,'[-]')
title('{\bfd.} pe-ratio 120m','fontweight','normal','fontsize',10)

nexttile(6)
% subplot(2,3,6)
 pe=nansum(Export_t(:,:,7,:),4)./nansum(nansum(NPP_t,3),4);
  pe(isinf(pe))=NaN;
h=function_plot_map_contours(pe,0:0.01:0.1);
caxis([0 0.1])
ylabel(h,'[-]')
title('{\bfe.} pe-ratio 1080m','fontweight','normal','fontsize',10)
colormap(viridis(100))
% colormap(flip(cmocean('Deep')))

% print -depsc figures_paper/fig3_rates_year



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPEPOD BIOMASS - Lopez anadon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dat=importdata('data_cop_biom_lopez/data_lopez_anadon.xlsx');
lat_dat=dat.data(:,10);
lon_dat=dat.data(:,9);
lon_dat=mod(lon_dat(:),360);
cop_biom_dat=dat.data(:,7);
naup_biom_dat=dat.data(:,8);
station_dat=dat.data(:,2);
cop_sfrac_dat=dat.data(:,1);
cop_body_mass=dat.data(:,5).*1e-3;
cop_mm=(cop_body_mass./exp(-16.41)).^(1/2.74)*1e-3;
cop_mm_mat=reshape(cop_mm,[length(cop_mm)/4,4]);
cop_biom_mat=reshape(cop_biom_dat,[length(cop_mm)/4,4]);
wght_mean=sum((cop_biom_mat.*cop_mm_mat),2)./sum(cop_biom_mat,2);

s1=find(cop_sfrac_dat==1); %Smallets
s2=find(cop_sfrac_dat==2);
s3=find(cop_sfrac_dat==3);
s4=find(cop_sfrac_dat==4); %Largest

sb4=find(param.Wvec>50);
sb3=find(param.Wvec>10 & param.Wvec<=50);
sb2=find(param.Wvec>1 & param.Wvec<=10);
sb1=find(param.Wvec<=1);

Cdtot_4D=matrixToGrid(sum(Cdeach,3), [], 'boxes.mat', 'grid.mat');
Cdtot_4D_1=matrixToGrid(sum(Cdeach(:,:,sb1),3), [], 'boxes.mat', 'grid.mat');
Cdtot_4D_2=matrixToGrid(sum(Cdeach(:,:,sb2),3), [], 'boxes.mat', 'grid.mat');
Cdtot_4D_3=matrixToGrid(sum(Cdeach(:,:,sb3),3), [], 'boxes.mat', 'grid.mat');
Cdtot_4D_4=matrixToGrid(sum(Cdeach(:,:,sb4),3), [], 'boxes.mat', 'grid.mat');

Cbiom_model=zeros(1,length(lat_dat(s1)));
Cbiom_model_1=zeros(1,length(lat_dat(s1)));
Cbiom_model_2=zeros(1,length(lat_dat(s1)));
Cbiom_model_3=zeros(1,length(lat_dat(s1)));
Cbiom_model_4=zeros(1,length(lat_dat(s1)));
for i=1:length(lat_dat(s1))
    idx1=find(y>lat_dat(i)-2.8 & y<lat_dat(i)+2.8);
    idx2=find(x>lon_dat(i)-3 & x<lon_dat(i)+3);
    Cbiom_model(i)=nansum(mean(Cdtot_4D(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
    Cbiom_model_1(i)=nansum(mean(Cdtot_4D_1(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
    Cbiom_model_2(i)=nansum(mean(Cdtot_4D_2(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
    Cbiom_model_3(i)=nansum(mean(Cdtot_4D_3(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
    Cbiom_model_4(i)=nansum(mean(Cdtot_4D_4(idx2(1),idx1(1),:,9*30+14:10*30+10).*permute(deltaz,[3,2,1]),4),3);
end

cop_tot_dat=cop_biom_dat(s1)+cop_biom_dat(s2)+cop_biom_dat(s3)+cop_biom_dat(s4);

%% Plot comparing data and model

size_mat=[cop_biom_dat(s4),cop_biom_dat(s3),cop_biom_dat(s2),cop_biom_dat(s1)+naup_biom_dat(s1)];
size_mat_mod=[Cbiom_model_4',Cbiom_model_3',Cbiom_model_2',Cbiom_model_1'];

x0=0;
y0=0;
width=18;
height=7;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(1,3,'Padding','Compact');

nexttile(2)
scatter(0,0,'markeredgecolor','none')
hold on
hb=bar(station_dat(s1),size_mat_mod./1000,'stacked');
cmap=flip(brewermap(4,'YlGnBu'));
hb(1).FaceColor = cmap(1,:); %large
hb(2).FaceColor = cmap(2,:);
hb(3).FaceColor = cmap(3,:);
hb(4).FaceColor = cmap(4,:); %small

leg1=legend([hb(4) hb(3) hb(2) hb(1)],{'<1','1 - 10','10 - 50','50<'});
legend boxoff
title(leg1,'Body-mass range')
set(gcf,'color','w');
ylim([0 10])
title('  {\bfa.} Model','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel('Biomass [gC m^{-2}]')
set(gca, 'XDir','reverse')
xlabel('Station')

nexttile(3)
scatter(0,0,'markeredgecolor','none')
hold on
hb=bar(station_dat(s1),size_mat./1000,'stacked');
cmap=flip(brewermap(4,'YlGnBu'));
hb(1).FaceColor = cmap(1,:);
hb(2).FaceColor = cmap(2,:);
hb(3).FaceColor = cmap(3,:);
hb(4).FaceColor = cmap(4,:);
hold on
hp=plot(station_dat(s1),(Cbiom_model)./1000,'k--','linewidth',1);
leg2=legend([hb(4) hb(3) hb(2) hb(1)],{'0.4','3.3','13.6','86.8'});
legend boxoff
title(leg2,'Average body-mass')
ylim([0 10])
title('  {\bfb.} Data','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Station')
set(gca, 'XDir','reverse')



c_p= [236 197 68]./255;
latlim = [-50 60];
lonlim = [-50 10];
nexttile(1)

axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on')
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
scatterm(lat_dat,lon_dat,10,'filled','markerfacecolor',c_p,'markeredgecolor','k')
text(lat_dat+2,lon_dat+2, num2str(1),'fontsize',8);
for i=1:24
    if i<10
        textm(lat_dat(i)+1,lon_dat(i)-5, num2str(i),'fontsize',6);
    else
        textm(lat_dat(i)+1,lon_dat(i)+2, num2str(i),'fontsize',6);        
    end

end

pos = get(gca, 'Position');
pos(1) = -0.45;
pos(3)=1.3;
set(gca, 'Position', pos)


set(gcf,'color','w');

% print -depsc figures_paper/fig3_copepods_data


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROTISTS BIOMASS 
%San Martin, E., Harris, R. P., & Irigoien, X. (2006). 
%Latitudinal variation in plankton size spectra in the Atlantic Ocean. 
%Deep Sea Research Part II: Topical Studies in Oceanography, 53(14-16), 1560-1572.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datP=importdata('data_protists/data_protists_biomass.xlsx');
lat_datP=datP.data(:,1);
lon_datP=datP.data(:,4);
lon_datP=mod(lon_datP(:),360);
month_datP=datP.data(:,6);
prot_biom_datP=datP.data(:,2);
prot_transect=datP.data(:,3);


Pdtot_4D=matrixToGrid(sum(Pdeach(:,:,param.V>2e-4),3), [], 'boxes.mat', 'grid.mat');
Pbiom_model=zeros(1,length(lat_datP));
spec_expo_model=zeros(1,length(lat_datP));
for i=1:length(lat_datP)
    
    idx1=find(y>lat_datP(i)-2.8 & y<lat_datP(i)+2.8);
    idx2=find(x>lon_datP(i)-3 & x<lon_datP(i)+3);
    Pbiom_model(i)=nansum(mean(Pdtot_4D(idx2(1),idx1(1),1,month_datP(i)*30-14-30:month_datP(i)*30+14-30).*deltaz(1),4),3);
    spec_expo_model(i)=(mean(Size_spec_expo(idx2(1),idx1(1),month_datP(i)*30-14-30:month_datP(i)*30+14-30),3));

end


%% plot protists biomass

cmap=[23, 25, 82;
64, 135, 163;
115, 209, 184]./255;

x0=0;
y0=0;
width=15;
height=9;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(3,3,'Padding','Compact');

nexttile(2,[1,2])
scatter(lat_datP(prot_transect==1),prot_biom_datP(prot_transect==1),15,'d','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==1),Pbiom_model(prot_transect==1)./1000,'color',cmap(1,:),'linewidth',1.5)
title('  {\bfa.} AMT 12 - June','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
set(gca,'box','off') 

nexttile(5,[1,2])
scatter(lat_datP(prot_transect==3),prot_biom_datP(prot_transect==3),15,'s','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==3),Pbiom_model(prot_transect==3)./1000,'color',cmap(2,:),'linewidth',1.5)
title('  {\bfb.} AMT 13 - October','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
set(gca,'box','off') 
ylabel('Biomass of nano and micro plankton [gC m^{-2}]')

nexttile(8,[1,2])
scatter(lat_datP(prot_transect==2),prot_biom_datP(prot_transect==2),15,'^','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==2),Pbiom_model(prot_transect==2)./1000,'color',cmap(3,:),'linewidth',1.5)
title('  {\bfc.} AMT 14 - May','fontweight','normal')
ax = gca;
xticks([-50:25:50])
ax.TitleHorizontalAlignment = 'left';
xlabel('Latitude')
set(gca,'box','off') 

set(gcf,'color','w');


nexttile(1,[3,1])
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on')
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
SC1=scatterm(lat_datP(prot_transect==1),lon_datP(prot_transect==1),15,'d','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC2=scatterm(lat_datP(prot_transect==3),lon_datP(prot_transect==3),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
SC3=scatterm(lat_datP(prot_transect==2),lon_datP(prot_transect==2),15,'^','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');

set(gcf,'color','w');
% print -depsc figures_paper/fig3_protists_data


%% plot size spectra data

datexpo=importdata('data_protists/data_exponent_lopezanadon.csv');
lat_dat_expo=datexpo.data(:,1);
expo_dat=datexpo.data(:,2);
amt_expo=datexpo.data(:,3);

cmap=[23, 25, 82;
64, 135, 163;
115, 209, 184]./255;

x0=0;
y0=0;
width=19;
height=5;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(1,3,'Padding','Compact');

nexttile
scatter(lat_dat_expo(amt_expo==12),expo_dat(amt_expo==12),15,'d','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==1),spec_expo_model(prot_transect==1),'color',cmap(1,:),'linewidth',1.5)
title('  {\bfa.} AMT 12 - June','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
ylabel('Exponent size spectrum [-]')
ylim([-1.3 -0.6])
xlabel('Latitude')

nexttile
scatter(lat_dat_expo(amt_expo==13),expo_dat(amt_expo==13),15,'s','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==3),spec_expo_model(prot_transect==3),'color',cmap(2,:),'linewidth',1.5)
title('  {\bfb.} AMT 13 - October','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
ylim([-1.3 -0.6])
xlabel('Latitude')

nexttile
scatter(lat_dat_expo(amt_expo==14),expo_dat(amt_expo==14),15,'^','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==2),spec_expo_model(prot_transect==2),'color',cmap(3,:),'linewidth',1.5)
title('  {\bfc.} AMT 14 - May','fontweight','normal')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
ylim([-1.3 -0.6])
xlabel('Latitude')

set(gcf,'color','w');

% print -depsc figures_paper/fig3_exponent_data


%% total figure prots and size spec

cmap=[23, 25, 82;
64, 135, 163;
115, 209, 184]./255;

msiz=8;

x0=0;
y0=0;
width=15;
height=8;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');
tiledlayout(3,3);%,'Padding','Compact');

nexttile(2)
scatter(lat_datP(prot_transect==1),prot_biom_datP(prot_transect==1),msiz,'d','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==1),Pbiom_model(prot_transect==1)./1000,'color',cmap(1,:),'linewidth',1)
text(0.05,0.95,'{\bfa.}','Units','normalized','fontsize',8)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
xticklabels([])


set(gca,'box','off') 

nexttile(5)
scatter(lat_datP(prot_transect==3),prot_biom_datP(prot_transect==3),msiz,'s','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==3),Pbiom_model(prot_transect==3)./1000,'color',cmap(2,:),'linewidth',1)
text(0.05,0.95,'{\bfc.}','Units','normalized','fontsize',8)
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xticks([-50:25:50])
xticklabels([])
set(gca,'box','off') 
ylabel({'Biomass of nano and';'micro plankton [gC m^{-2}]'})
ylim([0 12])

nexttile(8)
scatter(lat_datP(prot_transect==2),prot_biom_datP(prot_transect==2),msiz,'^','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==2),Pbiom_model(prot_transect==2)./1000,'color',cmap(3,:),'linewidth',1)
text(0.05,0.95,'{\bfe.}','Units','normalized','fontsize',8)
ax = gca;
xticks([-50:25:50])
ax.TitleHorizontalAlignment = 'left';
xlabel('Latitude')
set(gca,'box','off') 

set(gcf,'color','w');

nexttile(3)
scatter(lat_dat_expo(amt_expo==12),expo_dat(amt_expo==12),msiz,'d','markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==1),spec_expo_model(prot_transect==1),'color',cmap(1,:),'linewidth',1)
xticks([-50:25:50])
xticklabels([])
text(0.05,0.95,'{\bfb.}','Units','normalized','fontsize',8)
ylim([-1.3 -0.6])

nexttile(6)
scatter(lat_dat_expo(amt_expo==13),expo_dat(amt_expo==13),msiz,'s','markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==3),spec_expo_model(prot_transect==3),'color',cmap(2,:),'linewidth',1)
xticks([-50:25:50])
xticklabels([])
ylabel('Exponent size spectrum [-]')
text(0.05,0.95,'{\bfd.}','Units','normalized','fontsize',8)
ylim([-1.3 -0.6])

nexttile(9)
scatter(lat_dat_expo(amt_expo==14),expo_dat(amt_expo==14),msiz,'^','markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5)
hold on
plot(lat_datP(prot_transect==2),spec_expo_model(prot_transect==2),'color',cmap(3,:),'linewidth',1)
xticks([-50:25:50])
xlabel('Latitude')
text(0.05,0.95,'{\bff.}','Units','normalized','fontsize',8)
ylim([-1.3 -0.6])


nexttile(1,[3,1])
axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,...
    'Frame','on','Grid','on','MeridianLabel','on','ParallelLabel','on','fontsize',8)
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
axis off
SC1=scatterm(lat_datP(prot_transect==1),lon_datP(prot_transect==1),15,'d','filled','markerfacecolor',cmap(1,:),'markeredgecolor','none');
SC2=scatterm(lat_datP(prot_transect==3),lon_datP(prot_transect==3),15,'s','filled','markerfacecolor',cmap(2,:),'markeredgecolor','none');
SC3=scatterm(lat_datP(prot_transect==2),lon_datP(prot_transect==2),15,'^','filled','markerfacecolor',cmap(3,:),'markeredgecolor','none');
leg=legend([SC1,SC2,SC3],{'AMT 12','AMT 13','AMT 14'},'fontsize',8);
leg.ItemTokenSize(1) = 10;
leg.Location='southoutside';
leg.NumColumns=3;

% print -depsc figures_paper/fig3_exponent_data



%% export from fp and dead cells

Export_F_large_t=Export_F_t-Export_F_small_t;

x0=0;
y0=0;
width=18;
height=10;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');

nexttile
h=function_plot_map_contours(sum(Export_D_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.1:1);
ylabel(h,'[-]')
title({'{\bfa.} Export_D/Export_{tot} 120m'},'fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_small_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.1:1);
ylabel(h,'[-]')
title('{\bfb.} Export_{FP,small}/Export_{tot} 120m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_large_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.1:1);
ylabel(h,'[-]')
title('{\bfc.} Export_{FP,large}/Export_{tot} 120m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_D_t(:,:,8,:),4)./sum(Export_t(:,:,8,:),4),0:0.1:1);
ylabel(h,'[-]')
title('{\bfd.} Export_{D}/Export_{tot} 1080m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_small_t(:,:,8,:),4)./sum(Export_t(:,:,8,:),4),0:0.1:1);
ylabel(h,'[-]')
title('{\bfe.} Export_{FP,small}/Export_{tot} 1080m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours((sum(Export_F_large_t(:,:,8,:),4))./sum(Export_t(:,:,8,:),4),0:0.1:1);
ylabel(h,'[-]')
title('{\bff.} Export_{FP,large}/Export_{tot} 1080m','fontweight','normal','fontsize',10)
caxis([0 1])

colormap(viridis(100));

% print -depsc figures_paper/fig3_export_FP


%% yearly correlation
idx_cor=[1:8,10];

load('tempday.mat')
sst_t=matrixToGrid(tempday(:,1:2:end), [], 'boxes.mat', 'grid.mat');
lat_2D=repmat(y,[1,nx])';

z_exp_1=2;
z_exp_2=7;
ulat=70;%55;
lowlat=-60%-50;

bathy_surf=bathy(:,:,1);

cormat_yr=zeros(10,10,365);
rmse_ef=zeros(5,365);
totbiom=sum(Size_spec_sheldon,4);

Size_spec_coeff_year=mean(Size_spec_coeff,3);
Size_spec_expo_year=mean(Size_spec_expo,3);
Export_year_100=nansum(Export_t(:,:,z_exp_1,:),4);
Export_year_1000=nansum(Export_t(:,:,z_exp_2,:),4);
NPP_year=nansum(nansum(NPP_t,4),3);
pe_ratio_year_100=Export_year_100./nansum(nansum(NPP_t(:,:,1:2,:),4),3);
pe_ratio_year_1000=Export_year_1000./NPP_year;
sst_year=mean(sst_t(:,:,1,:),4);
totbiom_year=nanmean(totbiom,3);
Cbiom_year=nanmean(nansum(Cd_t,3),4);
TLyear=TLmean_year;%

TLyear(TLyear==0)=NaN;

Cbiom_year(isinf(Cbiom_year))=NaN;

Size_spec_coeff_year(bathy_surf==0)=NaN;
Size_spec_expo_year(bathy_surf==0)=NaN;
Export_year_100(bathy_surf==0)=NaN;
Export_year_1000(bathy_surf==0)=NaN;
NPP_year(bathy_surf==0)=NaN;
pe_ratio_year_100(bathy_surf==0)=NaN;
pe_ratio_year_1000(bathy_surf==0)=NaN;
sst_year(bathy_surf==0)=NaN;
totbiom_year(bathy_surf==0)=NaN;
Cbiom_year(bathy_surf==0)=NaN;

pe_ratio_year_100(pe_ratio_year_100>1)=NaN;
pe_ratio_year_1000(pe_ratio_year_1000>1)=NaN;

pe_ratio_year_100(pe_ratio_year_100==0)=NaN;
pe_ratio_year_1000(pe_ratio_year_1000==0)=NaN;

Export_year_100(Export_year_100==0)=NaN;
Export_year_1000(Export_year_1000==0)=NaN;

pe_ratio_year_100(isinf(pe_ratio_year_100))=NaN;
pe_ratio_year_1000(isinf(pe_ratio_year_1000))=NaN;

pe_ratio_year_100(isoutlier(pe_ratio_year_100))=NaN;
pe_ratio_year_1000(isoutlier(pe_ratio_year_1000))=NaN;

NPP_year(NPP_year==0)=NaN;

tbl = table(...
    log10(Size_spec_coeff_year(lat_2D>lowlat & lat_2D<ulat)),...
    Size_spec_expo_year(lat_2D>lowlat & lat_2D<ulat),...
    log10(NPP_year(lat_2D>lowlat & lat_2D<ulat)./1000),...
    log10(Export_year_100(lat_2D>lowlat & lat_2D<ulat)./1000),...
    log10(Export_year_1000(lat_2D>lowlat & lat_2D<ulat)./1000),...
    log10(pe_ratio_year_100(lat_2D>lowlat & lat_2D<ulat)),...
    log10(pe_ratio_year_1000(lat_2D>lowlat & lat_2D<ulat)),...
	sst_year(lat_2D>lowlat & lat_2D<ulat),...
    lat_2D(lat_2D>lowlat & lat_2D<ulat),...
    log10(Cbiom_year(lat_2D>lowlat & lat_2D<ulat)./1000),...
    TLyear(lat_2D>lowlat & lat_2D<ulat),...
    'VariableNames',{'int','slope','npp','exprt_1','exprt_2','ef_1','ef_2','tmp','lat','biom','TL'});

tbl_mat=table2array(tbl);


%% Final plot regressions

ylabs={{'Log_{10}Export_{120}';'[gC m^{-2} yr^{-1}]'},{'Log_{10}Export_{1080}';'[gC m^{-2} yr^{-1}]'},{'Log_{10}';'pe-ratio_{120} [-]'},{'Log_{10}';'pe-ratio_{1080} [-]'}};
xlabs={{'Exponent';'[-]'},{'Log_{10}Cop.biom.';'[gC m^{-2}]'},{'Trophic lvl.';'[-]'},{'Log_{10}NPP';'[gC m^{-2} yr^{-1}]'},{'SST';'[^{\circ}C]'}};

st=1;
x0=0;
y0=0;
width=18;
height=11;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

tiledlayout(4,5,'TileSpacing','Compact','Padding','Compact')

cmap=[23, 25, 82;
64, 135, 163;
135, 186, 201]./255;

strletters='abcdefghijklmnopqrst';

sty=1;
stx=1;
st=1;
for i=4:7  
    for j=[2,10,11,3,8]
        nexttile
        scatter(tbl_mat(:,j),tbl_mat(:,i),2,tbl_mat(:,9))
        lm=fitlm(tbl_mat(:,j),tbl_mat(:,i));
        hold on
%         plot(tbl_mat(:,j),lm.Coefficients.Estimate(1)+tbl_mat(:,j).*lm.Coefficients.Estimate(2),'k')
       
        if st<16
            xticklabels([])
        end
        if any([2:5:20, 3:5:20, 4:5:20, 5:5:20] == st)
            yticklabels([])
        end
        ylabel([])
        xlabel([])
        corval=corr([tbl_mat(:,j),tbl_mat(:,i)],'Rows','complete','type','Spearman');
        text(0.72,0.95,num2str(corval(1,2),'%.2f'),'Units','normalized','fontsize',8)
%         text(0.72,0.95,num2str(lm.Rsquared.Adjusted,'%.2f'),'Units','normalized','fontsize',8)
%         text(0.72,0.85,num2str(lm.RMSE,'%.4f'),'Units','normalized','fontsize',8)
%         text(0.62,0.95,strcat(['R^2=',num2str(lm.Rsquared.Adjusted,'%.2f')]),'Units','normalized','fontsize',8)
%         text(0.62,0.75,strcat(['RMSE=',num2str(lm.RMSE,'%.4f')]),'Units','normalized','fontsize',8)
        text(0.05,0.95,strcat(['\bf{',strletters(st),'}']),'Units','normalized','fontsize',8)
%         text(0.05,0.95,strcat(['\bf{',strletters(st),'}',...
%             '. R^2=',num2str(lm.Rsquared.Adjusted,'%.2f'),...
%             ', RMSE=',num2str(lm.RMSE,'%.4f')]),'Units','normalized','fontsize',8,'fontweight','normal')
        set(gca,'box','off') 
        colormap(brewermap(100,'RdBu'))
        caxis([-70 70])
        
        if j==2
            ylabel(ylabs{sty},'interpreter','tex')
            sty=sty+1;
        end
        
        if i==4
            ylim([0.5 2])
        elseif i==5
            ylim([-0.7 1.5])
        elseif i==6
            ylim([-1.5 -0])
%             yticks([0 0.25 0.5])
        elseif i==7
            ylim([-2.7 -0.8])
%             yticks([0 0.025 0.05])
        end
        %[2,10,11,3,8]
        if j==2
           xlim([-1.3 -0.6])      
        elseif j==10
           xlim([-1.3 1]) 
        elseif j==11
           xlim([2 4.2]) 
        elseif j==3
           xlim([1.5 2.6]) 
        elseif j==8
           xlim([-5 35]) 
        end
        
        if i==7
            xlabel(xlabs{stx},'interpreter','tex')
            stx=stx+1;
        end
                st=st+1;
    end
end
set(gcf,'color','w');
cb=colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Latitude')

% print -depsc figures_paper/fig3_corr


%% size spectrum slopes measured for only parts of the size spectrum

zs=1;%3 %idx of the chosen surface layer we integrate biomass 

%Pdeach is protists biomass with dimensions [nb,365,number of protists], where nb is the TM
%vector size 52749 (the same goes for Cdeach for copepods)
%deltaz_per is the depth bin size in m (I wanted depth integrated biomass of the surface layers)

Pd_surf=zeros(nx,ny,365,nbrP);
for i=1:nbrP
    Pspec1=matrixToGrid(Pdeach(:,:,i), [], 'boxes.mat', 'grid.mat');
    Pd_surf(:,:,:,i)=squeeze(nansum(Pspec1(:,:,1:zs,:).*deltaz_per(:,:,1:zs),3));    
end

Cd_surf=zeros(nx,ny,365,length(param.Wvec));
for i=1:length(param.Wvec)
    Cspec1=matrixToGrid(Cdeach(:,:,i), [], 'boxes.mat', 'grid.mat');
    Cd_surf(:,:,:,i)=squeeze(nansum(Cspec1(:,:,1:zs,:).*deltaz_per(:,:,1:zs),3));    
end

% Cd_surf(:)=NaN;

%make community size spectrum
bin_spec=logspace(-7,4,12*2); %boundaries of the community size spectrum
bin_spec_mean=geomean([bin_spec(1:end-1);bin_spec(2:end)]); %mean size for each bin
Size_spec=zeros(nx,ny,365,length(bin_spec_mean)); %normalized size spectrum
Size_specP=zeros(nx,ny,365,length(bin_spec_mean)); %normalized size spectrum
Size_specC=zeros(nx,ny,365,length(bin_spec_mean)); %normalized size spectrum
Size_spec_sheldon=zeros(nx,ny,365,length(bin_spec_mean)); %size spectrum with total biomass in each bin
delta_sizebin=bin_spec(2:end)-bin_spec(1:end-1);
% get total biomass in each size bin
for i=1:length(bin_spec_mean)
    binP=find(param.V>bin_spec(i) & param.V<bin_spec(i+1));
    binC=find(param.Wvec>bin_spec(i) & param.Wvec<bin_spec(i+1));
    Size_spec_sheldon(:,:,:,i)=(nansum(Cd_surf(:,:,:,binC),4)+nansum(Pd_surf(:,:,:,binP),4));
    Size_spec(:,:,:,i)=(nansum(Cd_surf(:,:,:,binC),4)+nansum(Pd_surf(:,:,:,binP),4))./delta_sizebin(i);
    Size_specP(:,:,:,i)=(nansum(Pd_surf(:,:,:,binP),4))./delta_sizebin(i);
    Size_specC(:,:,:,i)=(nansum(Cd_surf(:,:,:,binC),4))./delta_sizebin(i);
end

clearvars Pspec1 Cspec1 binP binC

%get parameters of the size spectrum (coeficient and exponent)
%takes some time to run...
Size_spec_expo=NaN(nx,ny,365);
Size_spec_coeff=NaN(nx,ny,365);

Size_spec_expo_P=NaN(nx,ny,365);
Size_spec_coeff_P=NaN(nx,ny,365);

Size_spec_expo_C=NaN(nx,ny,365);
Size_spec_coeff_C=NaN(nx,ny,365);

bin_spec_log=log10(bin_spec_mean)';
for ii=1:nx
    for  jj=1:ny
        for dd=1:365
            sp=Size_spec(ii,jj,dd,:);
            idx=isfinite(sp) & sp>0;
            if ~isempty(sp(idx))
                pp = polyfit(bin_spec_log(idx),log10(squeeze(sp(idx))),1);           
                Size_spec_expo(ii,jj,dd)=pp(1);
                Size_spec_coeff(ii,jj,dd)=10.^pp(2);
            end
            
            sp=Size_specP(ii,jj,dd,:);
            idx=isfinite(sp) & sp>0;
            if ~isempty(sp(idx))
                pp = polyfit(bin_spec_log(idx),log10(squeeze(sp(idx))),1);           
                Size_spec_expo_P(ii,jj,dd)=pp(1);
                Size_spec_coeff_P(ii,jj,dd)=10.^pp(2);
            end
            
            sp=Size_specC(ii,jj,dd,:);
            idx=isfinite(sp) & sp>0;
            if ~isempty(sp(idx))
                pp = polyfit(bin_spec_log(idx),log10(squeeze(sp(idx))),1);           
                Size_spec_expo_C(ii,jj,dd)=pp(1);
                Size_spec_coeff_C(ii,jj,dd)=10.^pp(2);
            end
        end
    end
end

Size_spec_expo(Size_spec_expo==0)=NaN;

figure
tiledlayout(1,3)
nexttile
h=function_plot_map_contours(mean(Size_spec_expo,3),-2:0.05:0);%,0:0.2:2.8)
ylabel(h,'\lambda [-]')
caxis([-1.5 -0.5])
% h.Ticks =-1.2:0.2:-0.6 ; 
title('{\bfa.} All sizes','fontweight','normal','fontsize',10)
colormap(brewermap(100,'BrBG'))

nexttile
h=function_plot_map_contours(mean(Size_spec_expo_P,3),-2:0.05:0);%,0:0.2:2.8)
ylabel(h,'\lambda [-]')
caxis([-1.5 -0.5])
% h.Ticks =-1.2:0.2:-0.6 ; 
title('{\bfb.} Only protists','fontweight','normal','fontsize',10)
colormap(brewermap(100,'BrBG'))

nexttile
h=function_plot_map_contours(mean(Size_spec_expo_C,3),-2:0.05:0);%,0:0.2:2.8)
ylabel(h,'\lambda [-]')
caxis([-1.5 -0.5])
% h.Ticks =-1.2:0.2:-0.6 ; 
title('{\bfc.} Only copepods','fontweight','normal','fontsize',10)
colormap(brewermap(100,'BrBG'))

% print -depsc figures_paper/fig3_spectrum_parts

%%

% figure; surface(Size_spec(:,:,1,1)'); shading flat

ii=1
jj=8
figure
tiledlayout(1,3)
for jij=1:3

%     if jij==1
%         ii=93;
%         jj=34;
%     elseif jij==2
%         ii=67;
%         jj=41;        
%     elseif jij==3
%         ii=63;
%         jj=48;        
%     end
    
    if jij==1
        ii=62;
        jj=51;
    elseif jij==2
        ii=119;
        jj=52;
    elseif jij==3
        ii=69;
        jj=40;
    end
    
    
    nexttile

bsl=10.^bin_spec_log;

for dd=1:1:365
sp=squeeze(Size_spec(ii,jj,dd,:));
plot(bsl,sp,'color',[214, 214, 214]./255)
hold on

end

sexpo=mean(Size_spec_expo(ii,jj,:),3);
scoeff=geomean(Size_spec_coeff(ii,jj,:),3);

sexpoP=mean(Size_spec_expo_P(ii,jj,:),3);
scoeffP=geomean(Size_spec_coeff_P(ii,jj,:),3);

sexpoC=mean(Size_spec_expo_C(ii,jj,:),3);
scoeffC=geomean(Size_spec_coeff_C(ii,jj,:),3);

sp=geomean(squeeze(Size_spec(ii,jj,:,:)),1);
h1=plot(bsl,(squeeze(sp)),'ko-','linewidth',1,'markersize',3,'markerfacecolor','k');
hold on
h2=plot(logspace(-7,4,50),scoeff.*logspace(-7,4,50).^sexpo,'r','linewidth',1.2);
h3=plot(logspace(-7,0,50),scoeffP.*logspace(-7,0,50).^sexpoP,'g','linewidth',1.2);
h4=plot(logspace(-1,4,50),scoeffC.*logspace(-1,4,50).^sexpoC,'b','linewidth',1.2);
set(gca,'yscale','log')
set(gca,'xscale','log')

if jij==1
    title('a. North Pacific')
    ylabel('mgC m^{-3} \mugC{-1}')
elseif jij==2
    title('b. North Atlantic')
elseif jij==3
    title('c. NP oligotrophic gyre')
    legend([h2,h3,h4],'All sizes','Only protists','Only copepods')
end

text(0.03,0.27,strcat(['\lambda = ',num2str(sexpo,'%.2f')]),'Units','normalized','fontsize',8)
text(0.03,0.17,strcat(['\lambda_P = ',num2str(sexpoP,'%.2f')]),'Units','normalized','fontsize',8)
text(0.03,0.07,strcat(['\lambda_C = ',num2str(sexpoC,'%.2f')]),'Units','normalized','fontsize',8)

xlabel('Body mass [mugC]')
ylim([1e-7 1e10])
xlim([1e-8 1e5])
yticks([1e-5 1e0 1e5 1e10])
xticks([1e-7 1e-3 1e1 1e5])

end

set(gcf,'color','w')

% print -depsc figures_paper/fig3_SI_spectrum_fit

%% fittin size spec for the seasonal plot

x1=119;
y1=52;

day_vec=[1 130 174 287];

ii=119;
jj=52;
figure
tiledlayout(2,2)
for jij=1:4

dd=day_vec(jij);
    nexttile

bsl=10.^bin_spec_log;


sp=squeeze(Size_spec(ii,jj,dd,:));
plot(bsl,sp,'color',[214, 214, 214]./255)
hold on



sexpo=Size_spec_expo(ii,jj,dd);
scoeff=Size_spec_coeff(ii,jj,dd);

sexpoP=Size_spec_expo_P(ii,jj,dd);
scoeffP=Size_spec_coeff_P(ii,jj,dd);

sexpoC=Size_spec_expo_C(ii,jj,dd);
scoeffC=Size_spec_coeff_C(ii,jj,dd);

sp=(squeeze(Size_spec(ii,jj,dd,:)));
h1=plot(bsl,(squeeze(sp)),'ko-','linewidth',1,'markersize',3,'markerfacecolor','k');
hold on
h2=plot(logspace(-7,4,50),scoeff.*logspace(-7,4,50).^sexpo,'r','linewidth',1.2);
h3=plot(logspace(-7,0,50),scoeffP.*logspace(-7,0,50).^sexpoP,'g','linewidth',1.2);
h4=plot(logspace(-1,4,50),scoeffC.*logspace(-1,4,50).^sexpoC,'b','linewidth',1.2);
set(gca,'yscale','log')
set(gca,'xscale','log')
if jij==1
    legend([h2,h3,h4],'All sizes','Only protists','Only copepods')
end

text(0.03,0.27,strcat(['\lambda = ',num2str(sexpo,'%.2f')]),'Units','normalized','fontsize',8)
text(0.03,0.17,strcat(['\lambda_P = ',num2str(sexpoP,'%.2f')]),'Units','normalized','fontsize',8)
text(0.03,0.07,strcat(['\lambda_C = ',num2str(sexpoC,'%.2f')]),'Units','normalized','fontsize',8)

ylim([1e-7 1e10])
xlim([1e-8 1e5])
yticks([1e-5 1e0 1e5 1e10])
xticks([1e-7 1e-3 1e1 1e5])

xlabel('Body mass [mugC]')
ylabel('Biomass spectrum\newline[mgC m^{-3} \mugC^{-1}]')
end

set(gcf,'color','w')

% print -depsc figures_paper/fig3_SI_spectrum_fit_seasonal

%%

figure
loglog(param.WD,param.sink_D,'k-o','linewidth',1,'MarkerFaceColor','k','markersize',4)
hold on
loglog(param.WF,param.sink_fp,'-s','color',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceColor',[0.5 0.5 0.5],'markersize',5)
xlabel('Particle carbon mass [\mugC]')
ylabel('Sinking rate [m d^{-1}]')
legend('Deadfalls','Fecal pellets')

set(gcf,'color','w')

%% locations of water columns etc

xid_fw=[63, 63, 63];
yid_fw=[40,43,49];

xid_season=[63, 119, 69];
yid_season=[49,52,40];

figure
scatter(wrapTo180(x(xid_season)),y(yid_season),35,'gd','filled','markeredgecolor','k');
hold on
scatter(wrapTo180(x(xid_fw)),y(yid_fw),40,'rp','filled','markeredgecolor','k');
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
set(gcf,'color','w');