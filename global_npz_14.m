%% Global model solved in water columns form- with several sizes of detritus
clearvars

% % 
c = parcluster('dcc R2019a');
c.AdditionalProperties.EmailAddress = 'mcsp@aqua.dtu.dk';
c.AdditionalProperties.MemUsage = '14GB';
c.AdditionalProperties.ProcsPerNode = 0;
c.AdditionalProperties.WallTime = '25:00'; 
c.saveProfile
% 
clust = parcluster('dcc R2019a');
numW=80;    % Exactly the number of nodes times the number of processors per cores requested
parpool(clust, numW);
% parpool('local',20);


yrs=20;%20; % years of the run

% Load Initial January TM and configuations
load('matrix_nocorrection_01.mat'); %loads transport matrix
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
load('grid.mat');
load('config_data.mat');
load('boxes.mat')

% Load paths for switching TM
loadPath = 'matrix_nocorrection_0';
loadPath1 =  'matrix_nocorrection_';

% %only needed if we change TM
% [idx_wc, bins]=function_find_idx_wc(Xbox, Ybox);
load('idx_wc.mat')
load('bins.mat')

% Preparing for timestepping. Using 43200 sec timesteps (0.5 days)
Ix = speye(nb,nb);
Aexp =Ix + (12*60*60)*Aexp;
Aimp = Aimp^(36);

Zbox_surf=Zbox;
zsurf=0;
for i=1:nz
   Zbox_surf(Zbox==z(i))=zsurf;
   zsurf=z(i);
end



%% Param setup
nbrP=7*2; %number of mixos size classes
nbrC_act=6; %number of copepods active feeders
nbrC_pass=3; %number of copepods passive feeders
nbrFP=8;
nbrD=8;

P_min=1e-7;
P_max=1e-1;
Cact_min=1.5;
Cact_max=10000;
Cpass_min=0.2;
Cpass_max=5;
C_size_classes=5;
totstatevars=1+nbrP+nbrC_act*C_size_classes+nbrC_pass*C_size_classes+nbrFP+nbrD;
param=parameters_copepod_model_5_clean(nbrP, P_min, P_max, nbrC_act, Cact_min, Cact_max,...
    nbrC_pass, Cpass_min, Cpass_max, C_size_classes, nbrFP, nbrD);

param.remin=0.04;%0.04; %remineralisation rate
param.sink_D=param.sink_D./5;%4;
param.pref_D(:)=1;
param.kw=0.03;


[theta_P_P,theta_P_D, theta_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=...
    func_feeding_kernels_2(param);

theta_cop_F(:)=1;
%% initial conditions

% % interpolated N data from world ocean atlas
N=load('N_oceandata.mat');
N=N.Nmati.*14;
N(:,:,8:end)=300;

% Uncomment if no spun up version is available
% if initial_conditions_switch==0
% N = zeros(nx,ny,nz);
P = zeros(nx,ny,nz);
C = zeros(nx,ny,nz);
D = zeros(nx,ny,nz);
F = zeros(nx,ny,nz);


% N(:,:,:) = 300;
P(:,:,:) = 0.01;
P(:,:,1:4) = 0.1;
C(:,:,:) = 0.01;
C(:,:,1:4) = 0.1;
% F(:,:,1)=0;

% Convert from grid form to vector form
N = gridToMatrix(N, [], 'boxes.mat', 'grid.mat');
P = gridToMatrix(P, [], 'boxes.mat', 'grid.mat');
C = gridToMatrix(C, [], 'boxes.mat', 'grid.mat');
D = gridToMatrix(D, [], 'boxes.mat', 'grid.mat');
F = gridToMatrix(F, [], 'boxes.mat', 'grid.mat');

Pmat=repmat(P,1,param.nbr_P);
Cmat=repmat(C,1,param.nbr_Ctot);
Fmat=repmat(F,1,param.nbr_fp);
Dmat=repmat(D,1,param.nbr_D);

month = 0;
imon=0;

% initialise vectors to save output (only for the last year)
years_save=1;
Nd = zeros(nb,365*years_save);
Pd = zeros(nb*nbrP,365*years_save);
Cd = zeros(nb*param.nbr_Ctot,365*years_save);
Dd = zeros(nb*param.nbr_D,365*years_save);
Fd = zeros(nb*param.nbr_fp,365*years_save);

% stuff
mon = [0 31 28 31 30 31 30 31 31 30 31 30 ];
mon2= [31 28 31 30 31 30 31 31 30 31 30 31];
ticks = (1+cumsum(mon))*2;
ticks(1) = 1;
labels = [{'J'} {'F'} {'M'} {'A'} {'M'} {'J'} {'J'} {'A'} {'S'} {'O'} {'N'} {'D'}];

xp = [x-x(1) ;360];

Tyear=365*2*yrs+1-365*2:2:365*2*yrs;
Tstart=365*2*yrs-365*2*years_save; 
%% Running the TM 
jj=1;
ti=1;
yrtime=1;

I_in=cell(1,length(bins));
mu_max_in=cell(1,length(bins));
AN_in=cell(1,length(bins));
remin_in=cell(1,length(bins));
AF_in=cell(1,length(bins));
AC_in=cell(1,length(bins));
R_in=cell(1,length(bins));
VN_in=cell(1,length(bins));
VF_in=cell(1,length(bins));
VL_in=cell(1,length(bins));

N_in=cell(1,length(bins));
Pmat_in3=cell(1,length(bins));
Cmat_in3=cell(1,length(bins));
Dmat_in=cell(1,length(bins));  
Fmat_in3=cell(1,length(bins));
Dprev2=cell(1,length(bins));
Fprev2=cell(1,length(bins));

N_out=cell(1,length(bins));
Pmat_out=cell(1,length(bins));
Cmat_out=cell(1,length(bins));
Fmat_out=cell(1,length(bins));
Dmat_out=cell(1,length(bins));

Ircell=cell(1,length(bins));

%idxB
idxBmort=cell(1,length(param.Wvec));
for i=1:length(param.Wvec)
    idxBmort{i}=find(param.Wvec>=param.Wvec(i)/10 & param.Wvec<=param.Wvec(i)*10);  
end

theta_opt=linspace(0,1,50)';

z0=[0; z];
deltaz=z0(2:end)-z0(1:end-1);
options1 = odeset('Refine',1);
diags=0;
refine=4;
options = odeset('Refine',refine);
opts1 = odeset('RelTol',1e-2,'AbsTol',1e-4);


load('parday.mat')
% idxday=0;
idxday=zeros(length(1:0.5:365*yrs),1);
ii=0;
for i=1:365*2*yrs
    ii=ii+1;
    idxday(i)=ceil(ii/2);    
    if mod(i,365*2)==0
        ii=0;
    end
end

tic

for i=1:365*2*yrs
    imon=imon+1;
    if imon>365*2
        imon=1;
    end
    % Test for time to change monthly TM
    if ismember(imon, ticks)
        month = month + 1;
        % Reset to Jan
        if month > 12
            month = 1;
            yrtime=yrtime+1;
            
        end
        % Load TM
        if month < 10  % el i es el doble de lo que deberia, m3 loads at 60
            load(strcat(loadPath, num2str(month), '.mat'));
            disp(strcat(loadPath, num2str(month), '.mat'));
        else
            load(strcat(loadPath1, num2str(month), '.mat'));
            disp(strcat(loadPath1, num2str(month), '.mat'));
        end
%         load(strcat('../../../woa_temperatures/T_data_month_',num2str(month),'.mat')) %load temp matrix
        load(strcat('T_',num2str(month),'.mat')); %load temp matrix
        T=gridToMatrix(Tmati, [], 'boxes.mat', 'grid.mat');
%         T(:)=10; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
            if ismember(imon, ticks)
               disp(['year=',num2str(yrtime), ' , month=', num2str(month)]);

            end
            
   
            
            
        Aexp=function_convert_TM_positive(Aexp);
        Aimp=function_convert_TM_positive(Aimp);
        
        % Preparing for timestepping. 43200s.
        Aexp =Ix + (12*60*60)*Aexp;
        Aimp = Aimp^(36);
          
     
    end
    
    N =  Aimp * ( Aexp  * N);
    % TM for each sizeclass of protist
    for ii=1:param.nbr_P        
        Pmat(:,ii) =  Aimp * ( Aexp  * Pmat(:,ii));       
    end
    
    for ii=1:param.nbr_Ctot        
        Cmat(:,ii) =  Aimp * ( Aexp  * Cmat(:,ii));       
    end
    
    for ii=1:param.nbr_fp        
        Fmat(:,ii) =  Aimp * ( Aexp  * Fmat(:,ii));       
    end
    
    for ii=1:param.nbr_D        
        Dmat(:,ii) =  Aimp * ( Aexp  * Dmat(:,ii));       
    end
    
    
                if any(isnan(N))
                  warning('there are NaNs just after the TM!!!!');
%                   break
                end 
    
    

    for j=1:length(bins)
        idx_1=idx_wc{j};

       T_inn=T(idx_1);

        %%%%%%%%%%%%%%%%%%%%%%%%%
       % temperature effects
       I_in{j}=param.I'.*param.Q_I.^((T_inn-param.Tref)/10);
       AN_in{j}=param.alpha_N.*param.Q_N.^((T_inn-param.Tref)/10);
       mu_max_in{j}=param.mu_max.*2.^((T_inn-18)/10);
       remin_in{j}=param.remin.*2.^((T_inn-param.Tref)/10);
       AF_in{j}=param.alpha_F.*1.5.^((T_inn-param.Tref)/10);
       AC_in{j}=param.alpha_c'.*1.5.^((T_inn-param.Tref)/10);
       R_in{j}=param.k'.*2.^((T_inn-param.Tref)/10);
       VN_in{j}=param.VN.*2.^((T_inn-18)/10);
       VF_in{j}=param.VF.*2.^((T_inn-18)/10);
       VL_in{j}=param.VL.*2.^((T_inn-18)/10);


        %%%%%%%%%%%%%%%%%%%%%%%%
        
%         N_in{idx_wcol(j,1),idx_wcol(j,2)}=permute(N(idx_wcol(j,1),idx_wcol(j,2),1:bins_wcol(idx_wcol(j,1),idx_wcol(j,2))),[3 2 1]);
        N_in{j}=N(idx_1);
        Pmat_in3{j}=Pmat(idx_1,:);
        Cmat_in3{j}=Cmat(idx_1,:);
        Dmat_in{j}=Dmat(idx_1,:);     
        Fmat_in3{j}=Fmat(idx_1,:);
% 
%         I0=par(idx_1);
%         P_layer=(sum(Pmat(idx_1,:),2).*deltaz(izBox(idx_1))); %(???)
%         ktot = kwdat(idx_1);%param.kw + param.kchl.*P_layer; % attenuation coeficient
%         Izh=I0.*exp(-ktot.*Zbox_surf(idx_1)); %light in the depth horizon
%         Il=Izh./(ktot.*deltaz(izBox(idx_1))).*(1-exp(-ktot.*deltaz(izBox(idx_1)))); %avg light in the mixed layer
%         Ircell{j}=Il;

        I0=parday(idx_1,idxday(i));
        P_layer=(sum(Pmat(idx_1,:),2).*deltaz(izBox(idx_1))); %(???)
        ktot = param.kw + param.kchl.*P_layer; % attenuation coeficient
        ktot_cumsum = param.kw + param.kchl.*cumsum(P_layer);
        Izh=I0.*exp(-ktot_cumsum.*Zbox_surf(idx_1)); %light starting in each depth horizon
        Il=Izh./(ktot.*deltaz(izBox(idx_1))).*(1-exp(-ktot.*deltaz(izBox(idx_1)))); %avg light experienced by phyto in each depth layer
        Ircell{j}=Il;
        

    end

%       ticBytes(gcp);  


%  ticBytes(gcp);
    parfor j=1:length(bins)
      idx_1=idx_wc{j};
      bins_wc=bins(j);
          options2 = odeset(options,opts1,'NonNegative',1:totstatevars*bins_wc);
          
          
          input_vars= [N_in{j}, Pmat_in3{j}, Cmat_in3{j}, Fmat_in3{j}, Dmat_in{j}];
          if any(isnan(input_vars(:)))
             [rownan, colnan]=find(isnan(input_vars))
              warning('there are NaNs in input_vars')
              input_vars(isnan(input_vars))=0;
              disp(j)
              
          end
          
   [t, Y] = ode45(@ode_mixo_17,[0 0.5], input_vars...
            , options2, param, theta_P_P, theta_cop_P, theta_cop_cop , theta_cop_F,...
               Ircell{j}, nbrP,I_in{j},AN_in{j},mu_max_in{j},...
              AF_in{j},AC_in{j},R_in{j},bins_wc,z,deltaz,remin_in{j},VN_in{j},VF_in{j},VL_in{j},...
              diags,theta_opt,idxBmort,Ybox(idx_1(1)),Xbox(idx_1(1)));

        if any(isnan(Y(:)))
           warning('there are nans after the odesolver')
           disp(j)
        end

        
        N_out{j} = Y(end,1:bins_wc);
        Pmat_out{j} = reshape(Y(end,bins_wc+1:bins_wc+bins_wc*param.nbr_P),[],nbrP);
        Cmat_out{j}= reshape(Y(end,bins_wc+bins_wc*param.nbr_P+1:bins_wc+bins_wc*param.nbr_P+bins_wc*param.nbr_Ctot),[],param.nbr_Ctot);
        Fmat_out{j} = reshape(Y(end,bins_wc+bins_wc*param.nbr_P+bins_wc*param.nbr_Ctot+1:end-bins_wc*param.nbr_D),[],param.nbr_fp);
        Dmat_out{j}= reshape(Y(end,end-bins_wc*param.nbr_D+1:end),[],param.nbr_D);
        
%         disp(j)

    end
% tocBytes(gcp);
      
    for j=1:length(bins)
      idx_1=idx_wc{j};
        N(idx_1) = N_out{j};
        Pmat(idx_1,:) = Pmat_out{j};
        Cmat(idx_1,:) = Cmat_out{j};
        Fmat(idx_1,:) = Fmat_out{j};
        Dmat(idx_1,:) = Dmat_out{j};
    end
     
   disp(i/2);
 
%start saving output for the last year run
    if i>=Tstart 
        if any(i==Tyear) % only saves every day (not every half day)
            Nd(:,ti) = N;
            Pd(:,ti) = Pmat(:);
            Cd(:,ti) = Cmat(:);
            Fd(:,ti) = Fmat(:);
            Dd(:,ti) = Dmat(:);
            ti=ti+1;
        end
    end



end
toc

%%

clearvars -except Nd Pd Cd Fd Dd t param nbrP layer yrs
% 
filename = 'ws_global_new_1.mat'; 
save(filename)
clearvars -except Pd Cd
save ws_Cd_new_1.mat Cd -v7.3
save ws_Pd_new_1.mat Pd -v7.3