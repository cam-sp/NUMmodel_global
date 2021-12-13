%% Stuff to load

clearvars
% 
load('ws_Pd_new_1_corr_z_12.mat')
load('ws_global_new_1_corr_z_12.mat')
load('ws_diagnostics_corr_z_12.mat')
load('ws_Cd_new_1_corr_z_12.mat')
load('frac_JLJR_diag_12.mat')


%%
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
%% create/clean diagnostic 4D matrices
Export = s_flux_D_diag + s_flux_F_diag;

NPP_year_integrated_1=nansum(nansum(NPP2_diag.*volb)).*1e-18;
Export_year_integrated_1=sum(sum(Export(izBox==2,:).*volb(izBox==2)./deltaz_vec(izBox==2))).*1e-18;
Export_year_integrated_2=sum(sum(Export(izBox==8,:).*volb(izBox==8)./deltaz_vec(izBox==8))).*1e-18;

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
% com_resp_t=matrixToGrid(com_resp_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
reproduction_t=matrixToGrid(reproduction_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
pred_C_Dtot_t=matrixToGrid(pred_C_Dtot_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
FPP_t=matrixToGrid(FPP_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
DPP_t=matrixToGrid(DPP_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
remin_t=matrixToGrid(remin_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;

NPPpico_t=matrixToGrid(NPP_pico_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
NPPnano_t=matrixToGrid(NPP_nano_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;
NPPmicro_t=matrixToGrid(NPP_micro_diag, [], 'boxes.mat', 'grid.mat').*deltaz_per;

frac_repr=squeeze(sum(reproduction_t(:,:,1:3,:),3)./sum(NPP_t(:,:,1:3,:),3));

frac_JLJR_t=matrixToGrid(frac_JLJR_diag, [], 'boxes.mat', 'grid.mat');


%% Biomasses

Pdeach=zeros(nb,365,nbrP);
st1=1;
st2=nb;
for i=1:nbrP
   Pdeach(:,:,i)=Pd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Pd_t=matrixToGrid(sum(Pdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total protists biomass


Cdeach=zeros(nb,365,param.nbr_Ctot);
st1=1;
st2=nb;
for i=1:param.nbr_Ctot
   Cdeach(:,:,i)=Cd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Cd_t=matrixToGrid(sum(Cdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total copepod biomass


% biomass detritus
Ddeach=zeros(nb,365,param.nbr_D);
st1=1;
st2=nb;
for i=1:param.nbr_D
   Ddeach(:,:,i)=Dd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Dd_t=matrixToGrid(sum(Ddeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total deadfalls biomass


Fdeach=zeros(nb,365,param.nbr_fp);
st1=1;
st2=nb;
for i=1:param.nbr_fp
   Fdeach(:,:,i)=Fd(st1:st2,:);
   st1=st1+nb;
   st2=st2+nb;
end
Fd_t=matrixToGrid(sum(Fdeach,3), [], 'boxes.mat', 'grid.mat').*deltaz_per; %total fecal pellets biomass


%% size spectrum bins

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
Size_spec_expo=zeros(nx,ny,365);
Size_spec_coeff=zeros(nx,ny,365);

bin_spec_log=log10(bin_spec_mean)';
for ii=1:nx
    for  jj=1:ny
        for dd=1:365
            pp = polyfit(bin_spec_log,log10(squeeze(Size_spec(ii,jj,dd,:))),1);           
            Size_spec_expo(ii,jj,dd)=pp(1);
            Size_spec_coeff(ii,jj,dd)=10.^pp(2);
        end
    end
end

% % %% length food-web 
% % 
% % % totbiom=sum(biom_bin.*permute(bin_spec(2:end)-bin_spec(1:end-1),[1,3,2]),3);
% % % cumsumbiom=cumsum(biom_bin.*permute(bin_spec(2:end)-bin_spec(1:end-1),[1,3,2]),3);
% % cumsum_sizespec=cumsum(Size_spec_sheldon,4);
% % totbiom=sum(Size_spec_sheldon,4);
% % biom_lengthFW=0.95*totbiom;
% % lengthFW=zeros(nx,ny,365);
% % for i=1:nx
% %   for  j=1:ny
% %       for id=1:365
% %         lengthFW_1=find(squeeze(cumsum_sizespec(i,j,id,:))>=biom_lengthFW(i,j,id));
% %         lengthFW(i,j,id)=bin_spec_mean(lengthFW_1(1));
% %       end
% %   end
% % end
% % lengthFW(squeeze(bathy(:,:,1))==0)=NaN;
% % 
% % %%
% % load tempday
% % %% average trophic level for all locations
% % 
% % TLmean=zeros(128,64);
% % TLmean_day=zeros(128,64,365);
% % 
% % for dayd=1:365
% % for xid=1:128
% %     for yid=1:64
% %         
% %         if bathy(xid,yid,2)==1
% % 
% % frac_auto=frac_JLJR_diag(izBox==1 & ixBox==xid & iyBox==yid,:,dayd);
% % Eflow_dep=zeros(54,54,2);
% % biomvec_dep=zeros(2,54);
% % for idep=1:2 %average over the 2 first depth ranges
% % 
% %     C=squeeze(Cdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
% %     P=squeeze(Pdeach(izBox==idep & ixBox==xid & iyBox==yid, dayd,:))'.*deltaz(idep);
% %     theta_P=theta_P_P;
% % 
% %     totsizevec=[param.V,param.Wvec'];
% % 
% %     % idx_prot=find(predsort_idx<=14);
% % 
% %     idx_parday=1:2:365*2;
% %     T_inn=tempday(izBox==1 & ixBox==xid & iyBox==yid,idx_parday(dayd));%T(idx_1);
% % 
% %             %%%%%%%%%%%%%%%%%%%%%%%%%
% %            % temperature effects
% %            I_in=param.I'.*param.Q_I.^((T_inn-param.Tref)/10);
% %            AN_in=param.alpha_N.*param.Q_N.^((T_inn-param.Tref)/10);
% %            mu_max_in=param.mu_max.*param.Q_k.^((T_inn-18)/10);
% %            AF_in=param.alpha_F.*1.5.^((T_inn-param.Tref)/10);
% %            alpha_c=param.alpha_c'.*1.5.^((T_inn-param.Tref)/10);
% %            VF_in=param.VF.*2.^((T_inn-18)/10);
% % 
% %     ingestion_each=zeros(length(param.Wvec),14+length(param.Wvec));
% %     contribution_preys_to_cops=zeros(length(param.Wvec),14+length(param.Wvec));
% %      for i=1:length(param.Wvec) %predators
% %         for j=1:14+length(param.Wvec) %preys           
% %             if j<=14 
% %                 E_C_tot=(theta_cop_P(i,j).*P(j));
% %             else
% %                 E_C_tot=(theta_cop_cop(i,j-14).*C(j-14));       
% %             end
% %             F_lvl=(alpha_c(i).*E_C_tot)./(alpha_c(i).*E_C_tot+I_in(i));
% %             ingestion_each(i,j)=I_in(i).*F_lvl.*C(i);             
% %         end
% %         ingestion_tot=sum(ingestion_each(i,:),2);
% %         for j=1:14+length(param.Wvec) %preys 
% %             contribution_preys_to_cops(i,j)=ingestion_each(i,j)./ingestion_tot;
% %         end
% %      end
% % 
% %     J_F=zeros(14,14);
% %     contribution_preys_to_prots=zeros(14,14);
% %      for i=1:14 %predators
% %         for j=1:i-1%size(E_P_tot,1) %preys            
% %             E_P_tot=P(j).*theta_P(i,j);
% %             J_F(i,j)=(VF_in(i).*AF_in(i).*E_P_tot)./(AF_in(i).*E_P_tot+VF_in(i)).*P(i);
% %         end
% %         ingestion_tot=sum(J_F(i,:),2);
% % 
% %         for j=1:i-1%size(E_P_tot,1) %preys               
% %             contribution_preys_to_prots(i,j)=J_F(i,j)./ingestion_tot;
% %         end
% % 
% %      end
% % 
% %     J_F2=[J_F,zeros(14,length(param.Wvec))];
% %     contribution_preys_to_prots_2=[contribution_preys_to_prots,zeros(14,length(param.Wvec))];
% %     contributions_dep(:,:,idep)=[contribution_preys_to_prots_2;contribution_preys_to_cops];
% %     Eflow_dep(:,:,idep)=[J_F2;ingestion_each];
% % 
% %     biomvec_dep(idep,:)=[P,C];
% % 
% % end
% % 
% % Eflow=mean(Eflow_dep,3); %flow of energy between compartments
% % biomvec=mean(biomvec_dep,1); %mean biomass in the two depth layers
% % contributions=mean(contributions_dep,3);
% % 
% % %to get the trophic level we need to organize all organisms relative to the
% % %sze of the prey they eat, since they all have different predator-prey mass
% % %ratios
% % preysize=[log10(param.V)-log10(param.beta_P),...
% %     log10(param.Wvec(param.ind_act)')-log10(param.beta_act),...
% %     log10(param.Wvec(param.ind_pass)')-log10(param.beta_pass)];
% % 
% % [predsort, predsort_idx]=sort(preysize); %now predators are sorted according to the size of the prey they eat
% % 
% % %fraction that is preffered to eat by each organisms
% % fracp=contributions'; % rows are preys and columns are predators
% % 
% % TL=ones(1,size(fracp,2));
% % TLpart=zeros(1,size(fracp,2));
% % 
% % frac_auto_2=[frac_auto,zeros(1,length(param.Wvec))];
% % 
% % %rows are prey, columns are predator
% % for i=predsort_idx%size(fracp,2) %predators
% %     for j=1:i-1%size(E_P_tot,1) %preys
% % 
% %                 TLpart(j)=TL(j).*fracp(j,i).*(1-frac_auto_2(i));
% % 
% %     end
% %     TL(i)=1+sum(TLpart);
% % end
% % 
% % TLmean(xid,yid)=sum(TL(14+param.ind_act).*biomvec(14+param.ind_act))./(sum(biomvec(14+param.ind_act)));
% % TLmean_day(xid,yid,dayd)=TLmean(xid,yid);
% % 
% %         end
% %     end
% % %     disp(xid)
% % end
% %     disp(dayd)
% % end
% % 
% % TLmean_year=nanmean(TLmean_day,3);

% % save('TLmean_year.mat','TLmean_year');

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
caxis([0 10])

ax1=nexttile;
h=function_plot_map_contours(nanmean(nansum(Cd_t(:,:,1:3,:),3),4)./1000,0:9);
ylabel(h,'[gC m^{-2}]')
title('{\bfb.} Copepods biomass','fontweight','normal','fontsize',10)
colormap(ax1,viridis(100));
caxis([0 10])

ax2=nexttile(4);
h=function_plot_map_contours(mean(Size_spec_expo,3),-1.4:0.05:-0.6);%,0:0.2:2.8)
ylabel(h,'[-]')
caxis([-1.3 -0.7])
h.Ticks =-1.2:0.2:-0.6 ; 
title('{\bfd.} Exponent (\lambda)','fontweight','normal','fontsize',10)
colormap(ax2, brewermap(100,'BrBG'))

ax3=nexttile(3);
h=function_plot_map_contours(TLmean_year,2:0.25:5.5);
caxis([2 4.5])
ylabel(h,'[-]')
title('{\bfc.} Trophic level active copepods','fontweight','normal','fontsize',10)
colormap(ax3,viridis(100));

print -depsc figures_paper/fig3_size_spec_global

%% biomass auto - mixo - hetero

frac_JLJR_diag_per=permute(frac_JLJR_diag,[1,3,2]);

frac_auto=zeros(size(frac_JLJR_diag_per));
frac_mixo=zeros(size(frac_JLJR_diag_per));
frac_hetero=zeros(size(frac_JLJR_diag_per));

autolim=0.7;
heterolim=0.3;

frac_auto(frac_JLJR_diag_per>=autolim)=1;
frac_mixo(frac_JLJR_diag_per>heterolim & frac_JLJR_diag_per<autolim)=1;
frac_hetero(frac_JLJR_diag_per<=heterolim)=1;

biom_auto_t=matrixToGrid(sum(Pdeach.*frac_auto,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_mixo_t=matrixToGrid(sum(Pdeach.*frac_mixo,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_hetero_t=matrixToGrid(sum(Pdeach.*frac_hetero,3), [], 'boxes.mat', 'grid.mat').*deltaz_per;

clearvars frac_JLJR_diag_per frac_auto frac_mixo frac_hetero

%% small cops biomass
biom_cops_small_t=matrixToGrid(sum(Cdeach(:,:,param.Wvec<=3),3), [], 'boxes.mat', 'grid.mat').*deltaz_per;
biom_cops_large_t=Cd_t-biom_cops_small_t;

%% seasoanl

Export_F_large_t=Export_F_t-Export_F_small_t;

ordervec=[1:3:18,2:3:18,3:3:18];
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

tiledlayout(7,3,'TileSpacing','Compact')%,'Padding','Compact');
% ha = tight_subplot(6,3,[.01 .05],[.1 .05],[.13 .12]);

for i=1:3
    if i==1
        x1=62;
        y1=51;
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
    legend boxoff
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
plot(pe_ratio_day_100,'-','color','k','linewidth',1.5)
hold on
plot(pe_ratio_day_1000,'-','color',[0.5 0.5 0.5],'linewidth',1.5)
xticklabels([])
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
set(gca,'box','off') 
xlim([0 365])
if i==3
    legend('120m','1000m')
    legend boxoff
end
if i==1
    ylabel({'[-]'})
end


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
ylim([-1.2 -0.6])
if i==1
    ylabel({'[-]'})
end

st=st+1;

end

set(gcf,'color','w');

print -depsc figures_paper/fig3_seasonal

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
height=12;
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
    legend([h1,h2,h3,h4,h5],{'Protists','Active cops.','Passive cops.','Community spec.\times100','OLS fit'})
    legend boxoff
end

end
set(gcf,'color','w');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dataset PP field data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars
NPP_int=squeeze(nansum(NPP_t,3));

NPP_int(NPP_int<1)=NaN;

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

figure
subplot(2,2,1)
scatter(NPP_tmm(NPP_tmm_std<0.1.*NPP_tmm),NPP_int(NPP_tmm_std<0.1.*NPP_tmm),10,lat_3D(NPP_tmm_std<0.1.*NPP_tmm));
hold on
plot(linspace(0,1500),linspace(0,1500),'k--')
colorbar
xlabel('Data [mgC m^{-2} day^{-1}]')
ylabel('Model [mgC m^{-2} day^{-1}]')
title('NPP')

subplot(2,2,2)
histogram(log10(NPP_tmm(NPP_tmm_std<0.1.*NPP_tmm)));
hold on
histogram(log10(NPP_int(~isnan(NPP_tmm) & NPP_tmm_std<0.1.*NPP_tmm)));
legend('Data','Model')
xlabel('Log_{10} NPP [mgC m^{-2} day^{-1}]')
ylabel('[]')
title('Log_{10} NPP')

subplot(2,1,2)
scatter(lon_3D(NPP_tmm_std<0.1.*NPP_tmm),lat_3D(NPP_tmm_std<0.1.*NPP_tmm),8,'g','filled');
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
title('NPP data locations')
set(gcf,'color','w');

idx_SO=find(lat_3D<-50);
idx_HOT=find(lon_3D<-150 & lat_3D>0);
idx_BATS=find(lon_3D>-67 & lon_3D<-65 & lat_3D>29 & lat_3D<30);
idx_NA=find(lon_3D>-21 & lon_3D<0);
idx_arab=find(lon_3D>50 & lon_3D<67);


NPP_int_log=log10(NPP_int(:));
NPP_tmm_log=log10(NPP_tmm(:));
NPP_int_log(isinf(NPP_int_log))=NaN;
NPP_tmm_log(isinf(NPP_int_log))=NaN;
NPP_int_log(isinf(NPP_tmm_log))=NaN;
NPP_tmm_log(isinf(NPP_tmm_log))=NaN;

NPPso_nan=NPP_tmm(idx_SO);
NPPso_nan=NPPso_nan(~isnan(NPPso_nan));
NPPhot_nan=NPP_tmm(idx_HOT);
NPPhot_nan=NPPhot_nan(~isnan(NPPhot_nan));
NPPbats_nan=NPP_tmm(idx_BATS);
NPPbats_nan=NPPbats_nan(~isnan(NPPbats_nan));
NPPna_nan=NPP_tmm(idx_NA);
NPPna_nan=NPPna_nan(~isnan(NPPna_nan));
NPPar_nan=NPP_tmm(idx_arab);
NPPar_nan=NPPar_nan(~isnan(NPPar_nan));

RMSD1=(1./length(NPPso_nan).*(nansum((log10(NPP_int(idx_SO))-log10(NPP_tmm(idx_SO))).^2))).^0.5;
RMSD2=(1./length(NPPhot_nan).*(nansum((log10(NPP_int(idx_HOT))-log10(NPP_tmm(idx_HOT))).^2))).^0.5;
RMSD3=(1./length(NPPbats_nan).*(nansum((log10(NPP_int(idx_BATS))-log10(NPP_tmm(idx_BATS))).^2))).^0.5;
RMSD4=(1./length(NPPna_nan).*(nansum((log10(NPP_int(idx_NA))-log10(NPP_tmm(idx_NA))).^2))).^0.5;
RMSD5=(1./length(NPPar_nan).*(nansum((log10(NPP_int(idx_arab))-log10(NPP_tmm(idx_arab))).^2))).^0.5;
RMSD=[RMSD1;RMSD2;RMSD3;RMSD4;RMSD5];


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

lon2_flux=mod(lon_data_flux(:),360);

zface=cumsum(deltaz);

mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
monticks = [(1+cumsum(mon)) 365];

NPP_t_100=squeeze(nansum(NPP_t(:,:,1:2,:),3));

flux_model=zeros(length(lon_data_flux),1);
eff_model=zeros(length(lon_data_flux),1);
PP_model_flux=zeros(length(lon_data_flux),1);
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
        eff_model(i)=flux_model(i)./nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
        PP_model_flux(i)=nanmean(NPP_t_100(xid(1),yid(1),day_data_flux(i)),3);
    end

    
end


idx_data=find(~isnan(eff_model) & ~isinf(eff_model) & ~isnan(eratio_data) & PP_model_flux>=0 & flux_data>0);

flux_model2=flux_model(idx_data);
eff_model2=eff_model(idx_data);
PP_model_flux2=PP_model_flux(idx_data);
flux_data2=flux_data(idx_data);
eratio2=eratio_data(idx_data);
lat_data2=lat_data_flux(idx_data);
lon_data2=lon_data_flux(idx_data);
PP_data_flux2=PP_data_flux(idx_data);


figure
subplot(1,2,1)
scatter(flux_data2,flux_model2,10,lat_data2)
hold on
plot(linspace(0,800),linspace(0,800),'k--')
xlabel('data')
ylabel('model')
set(gca,'yscale','log')
set(gca,'xscale','log')
colorbar
title('Carbon export')
ylim([1e1 1e3])
xlim([1e1 1e3])

subplot(1,2,2)
histogram(log10(flux_data2))
hold on
histogram(log10(flux_model2))
legend('data','model')
set(gcf,'color','w');

corr(log10(flux_data2),log10(flux_model2),'type','Spearman')


%% binned carbon export

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
scatter(lon_data2,lat_data2,8,'r','filled');
hold on
geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'black');
title('Carbon export data location')
set(gcf,'color','w');


%% with binned data

st=1;
x0=0;
y0=0;
width=18;
height=7;
fig=figure('Units','centimeters',...
'Position',[x0 y0 width height],...
'PaperPositionMode','auto');

greycolor=[0.5 0.5 0.5]

subplot(1,2,1)
ssiz=10;
cmap=brewermap(6,'Dark2');
errorbar(NPP_tmm(idx_SO),NPP_int(idx_SO),NPP_tmm_std(idx_SO),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
hold on
errorbar(NPP_tmm(idx_HOT),NPP_int(idx_HOT),NPP_tmm_std(idx_HOT),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
errorbar(NPP_tmm(idx_BATS),NPP_int(idx_BATS),NPP_tmm_std(idx_BATS),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
errorbar(NPP_tmm(idx_NA),NPP_int(idx_NA),NPP_tmm_std(idx_NA),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
errorbar(NPP_tmm(idx_arab),NPP_int(idx_arab),NPP_tmm_std(idx_arab),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
h1=scatter(NPP_tmm(idx_SO),NPP_int(idx_SO),ssiz,'markerfacecolor',cmap(3,:),'markeredgecolor',cmap(3,:),'markerfacealpha',0.5);
h2=scatter(NPP_tmm(idx_HOT),NPP_int(idx_HOT),ssiz,'markerfacecolor',cmap(2,:),'markeredgecolor',cmap(2,:),'markerfacealpha',0.5);
h3=scatter(NPP_tmm(idx_BATS),NPP_int(idx_BATS),ssiz,'markerfacecolor',cmap(1,:),'markeredgecolor',cmap(1,:),'markerfacealpha',0.5);
h4=scatter(NPP_tmm(idx_NA),NPP_int(idx_NA),ssiz,'markerfacecolor',cmap(4,:),'markeredgecolor',cmap(4,:),'markerfacealpha',0.5);
h5=scatter(NPP_tmm(idx_arab),NPP_int(idx_arab),ssiz,'markerfacecolor',cmap(5,:),'markeredgecolor',cmap(5,:),'markerfacealpha',0.5);
plot(linspace(0,2000),linspace(0,2000),'--','Color',[0.5 0.5 0.5])
legend([h1,h2,h3,h4,h5],'SO','HOT','BATS','NA','Arab')
set(gcf,'color','w');
xlabel('Data [mgC m^{-2} day^{-1}]','fontsize',10)
ylabel('Model [mgC m^{-2} day^{-1}]','fontsize',10)
title('{\bfa.} NPP by regions','fontweight','normal')
corval=corr(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)),'type','Spearman','Rows','complete');
% text(0.05,0.95,strcat(['r=',num2str(corval,'%.2f')]),'Units','normalized','fontsize',9)
rmse= fitlm(NPP_tmm(~isnan(NPP_tmm)),NPP_int(~isnan(NPP_tmm)));
Err= (1/length(~isnan(NPP_tmm)))*(nansum((NPP_int(:))-(NPP_tmm(:))));
ax = gca;
% ax.TitleHorizontalAlignment = 'left';
set(gca,'box','off') 

subplot(1,2,2)
errorbar(flux_data_tmm(flux_data_tmm_std~=0),flux_model(flux_data_tmm_std~=0),flux_data_tmm_std(flux_data_tmm_std~=0),'horizontal','o'...
    ,'markerfacecolor',greycolor,'markeredgecolor',greycolor,'color',greycolor,'markersize',1,'CapSize',2)
hold on
scatter(flux_data_tmm(:),flux_model(:),10,sst_3d(:),'filled','markerfacealpha',0.5)
plot(logspace(0,3),logspace(0,3),'--','color',[0.5 0.5 0.5])
colormap(viridis(100))
ch=colorbar;
ylabel(ch,'Latitude')
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([10^0.5 1e3])
xlim([10^0.5 1e3])
set(gcf,'color','w');
xlabel('Data [mgC m^{-2} day^{-2}]')
ylabel('Model [mgC m^{-2} day^{-2}]')
title('{\bfb.} Carbon export','fontweight','normal')
ax = gca;
set(gca,'box','off') 


print -depsc figures_paper/fig3_NPP_data


%% rates compared to data plots


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
h=function_plot_map_contours(nansum(nansum(NPP_t,3),4)./1000,0:20:500);
ylabel(h,'[gC m^{-2} year^{-1}]')
% caxis([0 9])
title('{\bfa.} NPP','fontweight','normal','fontsize',10)
caxis([0 500])

nexttile
% subplot(2,3,2)
h=function_plot_map_contours(nansum(Export_t(:,:,2,:),4)./1000,0:2:50);
ylabel(h,'[gC m^{-2} year^{-1}]')
title('{\bfb.} Export 120m','fontweight','normal','fontsize',10)
caxis([0 50])

nexttile
% subplot(2,3,3)
h=function_plot_map_contours(nansum(Export_t(:,:,7,:),4)./1000,0:2:50);
ylabel(h,'[gC m^{-2} year^{-1}]')
title('{\bfc.} Export 1080m','fontweight','normal','fontsize',10)
caxis([0 50])

nexttile(5)
% subplot(2,3,5)
 pe=nansum(Export_t(:,:,2,:),4)./nansum(nansum(NPP_t,3),4);
  pe(isinf(pe))=NaN;
h=function_plot_map_contours(pe,0:0.02:0.4);
caxis([0 0.4])
ylabel(h,'[-]')
title('{\bfd.} pe-ratio 120m','fontweight','normal','fontsize',10)

nexttile(6)
% subplot(2,3,6)
 pe=nansum(Export_t(:,:,7,:),4)./nansum(nansum(NPP_t,3),4);
  pe(isinf(pe))=NaN;
h=function_plot_map_contours(pe,0:0.01:0.2);
caxis([0 0.2])
ylabel(h,'[-]')
title('{\bfe.} pe-ratio 1080m','fontweight','normal','fontsize',10)
colormap(viridis(100))

print -depsc figures_paper/fig3_rates_year

% %% Plots NPP and export 
% 
% load('tempday.mat')
% sst_t=matrixToGrid(tempday(:,1:2:end), [], 'boxes.mat', 'grid.mat');
% lat_2D=repmat(y,[1,nx])';
% 
% x0=0;
% y0=0;
% width=17;
% height=11;
% fig=figure('Units','centimeters',...
% 'Position',[x0 y0 width height],...
% 'PaperPositionMode','auto');
% tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
% 
% nexttile
% h=function_plot_map_contours(mean(nansum(NPP_t,3),4));
% ylabel(h,'[mgC m^{-2} day^{-1}]')
% % caxis([0 9])
% title('{\bfa.} NPP','fontweight','normal','fontsize',10)
% caxis([0 1200])
% 
% nexttile
% h=function_plot_map_contours(mean(Export_t(:,:,2,:),4));
% % caxis([0 9])
% ylabel(h,'[mgC m^{-2} day^{-1}]')
% title('{\bfb.} Export 120m','fontweight','normal','fontsize',10)
% 
% 
% nexttile
% h=function_plot_map_contours(mean(Export_t(:,:,8,:),4));
% % caxis([-1 3.5])
% % h.Ticks =-1:3 ; 
% % h.TickLabels = {'10^{-1}','10^{0}','10^1','10^{2}','10^3'} ; 
% ylabel(h,'[mgC m^{-2} day^{-1}]')
% title('{\bfc.} Export 1080m','fontweight','normal','fontsize',10)
% 
% nexttile
% h=function_plot_map_contours(mean(sst_t(:,:,1,:),4),-5:2:40);
% ylabel(h,'[^{\circ}C]')
% % caxis([0 9])
% title('{\bfd.} SST','fontweight','normal','fontsize',10)
% 
% nexttile
% pe2y=nanmean(Export_t(:,:,2,:)./nansum(NPP_t(:,:,1:2,:),3),4);
% pe2y(isinf(pe2y))=NaN;
% h=function_plot_map_contours(pe2y,0:0.02:0.35);
% caxis([0 0.35])
% ylabel(h,'[-]')
% title('{\bfe.} pe-ratio 120m','fontweight','normal','fontsize',10)
% 
% 
% nexttile
% pe7y=nanmean(Export_t(:,:,7,:)./nansum(NPP_t,3),4);
% pe7y(isinf(pe7y))=NaN;
% h=function_plot_map_contours(pe7y,0:0.005:0.1);
% caxis([0 0.1])
% % h.Ticks =-1:3 ; 
% % h.TickLabels = {'10^{-1}','10^{0}','10^1','10^{2}','10^3'} ; 
% ylabel(h,'[-]')
% title('{\bff.} pe-ratio 1080m','fontweight','normal','fontsize',10)
% 
% 
% colormap(viridis(100))
% % colormap(flip(brewermap(100,'YlGnBu')))

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
height=8;
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

legend([hb(4) hb(3) hb(2) hb(1)],{'<1','1 - 10','10 - 50','50<'})
legend boxoff
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
legend([hb(4) hb(3) hb(2) hb(1)],{'0.4','3.3','13.6','86.8'})
legend boxoff
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
scatterm(lat_dat,lon_dat,15,'filled','markerfacecolor',c_p,'markeredgecolor','k')
text(lat_dat+2,lon_dat+2, num2str(1),'fontsize',8);
for i=1:24
    if i<10
        textm(lat_dat(i)+1,lon_dat(i)-5, num2str(i),'fontsize',8);
    else
        textm(lat_dat(i)+1,lon_dat(i)+2, num2str(i),'fontsize',8);        
    end

end

pos = get(gca, 'Position');
pos(1) = -0.45;
pos(3)=1.3;
set(gca, 'Position', pos)

set(gcf,'color','w');
print -depsc figures_paper/fig3_copepods_data


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
print -depsc figures_paper/fig3_protists_data


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
tiledlayout(3,3,'Padding','Compact');

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
h=function_plot_map_contours(sum(Export_D_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.05:1);
ylabel(h,'[-]')
title({'{\bfa.} Export_D/Export_{tot} 120m'},'fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_small_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.05:1);
ylabel(h,'[-]')
title('{\bfb.} Export_{FP,small}/Export_{tot} 120m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_large_t(:,:,2,:),4)./sum(Export_t(:,:,2,:),4),0:0.05:1);
ylabel(h,'[-]')
title('{\bfc.} Export_{FP,large}/Export_{tot} 120m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_D_t(:,:,8,:),4)./sum(Export_t(:,:,8,:),4),0:0.05:1);
ylabel(h,'[-]')
title('{\bfd.} Export_{D}/Export_{tot} 1080m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours(sum(Export_F_small_t(:,:,8,:),4)./sum(Export_t(:,:,8,:),4),0:0.05:1);
ylabel(h,'[-]')
title('{\bfe.} Export_{FP,small}/Export_{tot} 1080m','fontweight','normal','fontsize',10)
caxis([0 1])

nexttile
h=function_plot_map_contours((sum(Export_F_large_t(:,:,8,:),4))./sum(Export_t(:,:,8,:),4),0:0.05:1);
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
    (pe_ratio_year_100(lat_2D>lowlat & lat_2D<ulat)),...
    (pe_ratio_year_1000(lat_2D>lowlat & lat_2D<ulat)),...
	sst_year(lat_2D>lowlat & lat_2D<ulat),...
    lat_2D(lat_2D>lowlat & lat_2D<ulat),...
    log10(Cbiom_year(lat_2D>lowlat & lat_2D<ulat)./1000),...
    TLyear(lat_2D>lowlat & lat_2D<ulat),...
    'VariableNames',{'int','slope','npp','exprt_1','exprt_2','ef_1','ef_2','tmp','lat','biom','TL'});

tbl_mat=table2array(tbl);


%% good final plot regression

ylabs={{'Log_{10}Export_{120}';'[gC m^{-2} yr^{-1}]'},{'Log_{10}Export_{1080}';'[gC m^{-2} yr^{-1}]'},{'pe-ratio_{120}';'[-]'},{'pe-ratio_{1080}';'[-]'}};
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
        plot(tbl_mat(:,j),lm.Coefficients.Estimate(1)+tbl_mat(:,j).*lm.Coefficients.Estimate(2),'k')
       
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
        text(0.05,0.95,strcat(['\bf{',strletters(st),'}']),'Units','normalized','fontsize',8)
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
            ylim([0 1.5])
        elseif i==6
            ylim([0 0.5])
        elseif i==7
            ylim([0 0.15])
        end
        %[2,10,11,3,8]
        if j==2
           xlim([-1.2 -0.6])      
        elseif j==10
           xlim([-1 1]) 
        elseif j==11
           xlim([2 4.2]) 
        elseif j==3
           xlim([1.5 2.5]) 
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



        
