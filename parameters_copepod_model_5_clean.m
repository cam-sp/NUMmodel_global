function param=parameters_copepod_model_5_clean(P_size_classes, P_min, P_max, Cact_populations, Cact_min, Cact_max,...
    Cpass_populations, Cpass_min, Cpass_max, C_size_classes, nbrFp, nbrD)

% Here we create the whole set of parameters needed to run the model
%note that te names of parameters here are different than in the paper or
%the ones you introduced initially. But everything has been commented. If
%you have any questions feel free to email me: mcsp@aqua.dtu.dk

%Camila Serra 12/03/2020
%-------------------------------------------------------------------

%Make grid of protists
param.nbr_P=P_size_classes; %number of size-classes
V_points=param.nbr_P+1; %number of points in the grid
param.min_P=log10(P_min); %minimum size generalists in log10
param.max_P=log10(P_max); %maximum size generalists in log10
param.V_grid=logspace(param.min_P,param.max_P,V_points); %generalists grid
param.V_up=param.V_grid(2:end); %upper sizes
param.V_dw=param.V_grid(1:end-1); %lower sizes
param.V=geomean([param.V_up;param.V_dw]); %central sizes
param.delta_V=(param.V_up)-(param.V_dw); %bin size
param.ratio_V=(param.V_up./param.V_dw); %bin size

%Make grid of copepods
C_sp_act=Cact_populations; %number of species (populations) active feeders
C_sp_pass=Cpass_populations;%5; %number of species (population) passive feeders
param.C_sp_act=C_sp_act;
param.C_sp_pass=C_sp_pass;
param.C_sp=C_sp_act+C_sp_pass; %total number of populations
param.nbr_discr=C_size_classes;
param.min_cop_act=log10(Cact_min); % minimum size of adult active copepods in log10
param.max_cop_act=log10(Cact_max); % maximum size of adult active copepods in log10
param.min_cop_pass=log10(Cpass_min); % minimum size of adult passive copepods in log10
param.max_cop_pass=log10(Cpass_max); % maximum size of adult passive copepods in log10
C_grid_act=zeros(param.nbr_discr,C_sp_act); % each column is a population of active copepods, 
                                            % each row is a size bin, from juveniles (top) to adults (bottom)
C_grid_pass=zeros(param.nbr_discr,C_sp_pass); %same as above for passive copepods
Ca_sp_act=logspace(param.min_cop_act,param.max_cop_act,C_sp_act); %vector of adult copepod size
Ca_sp_pass=logspace(param.min_cop_pass,param.max_cop_pass,C_sp_pass); %vector of adult copepod size
OA_ratio_act=0.01;%0.002; %offspring to adult ratio 0.003 0.0063
OA_ratio_pass=0.01;%0.02; %offspring to adult ratio 0.01
for i=1:C_sp_act
    C_grid_act(:,i)=logspace(log10(Ca_sp_act(i).*OA_ratio_act),log10(Ca_sp_act(i)),param.nbr_discr); %copepod grid 0.003
end
for i=1:C_sp_pass
    C_grid_pass(:,i)=logspace(log10(Ca_sp_pass(i).*OA_ratio_pass),log10(Ca_sp_pass(i)),param.nbr_discr); %copepod grid 0.003
end

C_grid=cat(2,C_grid_act,C_grid_pass);
C_up=C_grid(2:end,:); %upper sizes
C_dw=C_grid(1:end-1,:); %lower sizes
% W=(C_up-C_dw)./log(C_up./C_dw); %Central sizes
W=sqrt(C_up.*C_dw); %Central sizes (geomean);
param.W=cat(1,W,C_up(end,:)); % total grid of copepod weights, last row is the weight of the adult
param.C_up=cat(1,C_up,C_up(end,:)); %we fix it to include the aduts, useful for the feeding kernels afterwards
param.C_dw=cat(1,C_dw,C_up(end,:)); %we fix it to include the aduts, useful for the feeding kernels afterwards
param.Wa=param.W(end,:); %adults mass
param.nbr_cops=length(param.Wa);
param.deltaC=param.C_up-param.C_dw;
param.nbr_stages=length(param.W(:,1));
% z=C_grid(1:end-1,:)./C_grid(2:end,:);
z=param.C_dw./param.C_up;
param.z=z(:);
param.nbr_Ctot=length(param.W(:));

param.ind_a=param.nbr_stages:param.nbr_stages:param.nbr_stages*param.nbr_cops; %index of adult copepods
param.ind_b=1:param.nbr_stages:param.nbr_stages*param.nbr_cops; %index of newborn nauplii
indexs=1:param.nbr_stages*param.nbr_cops;
param.ind_j=setdiff(indexs,param.ind_a); %indexs of all juveniles stages
param.ind_rest=setdiff(param.ind_j,param.ind_b); %indexs of juveniles stages except newborns
param.ind_act=1:param.nbr_stages*C_sp_act; %indexs of active copepods
param.ind_pass=param.ind_act(end)+1:length(param.W(:)); %index of passive copepods
param.Wvec=param.W(:); % convert the matrix into a column vector

%%%%%%%%%%%%%%%%%%%%%%%%%

%grid for the Fecal pellets size bins
FPV=10^4.547.*param.Wvec.^1; %find volume fecal pellets
mfp=10^(-0.665)*FPV.^(0.939) *1e-6; %mass of fecal pellets

param.nbr_fp=nbrFp; %number of bins
fp_vec=logspace(log10(min(mfp))-0.1*min(mfp), log10(max(mfp)+0.1*max(mfp)) ,param.nbr_fp+1);
max_len_fp=zeros(1,length(fp_vec)-1);
for  i=1:param.nbr_fp
    max_len_fp(i)=length(find(mfp>=fp_vec(i)  & mfp<fp_vec(i+1))); 
end

vec_idx_fp=zeros(1,length(param.Wvec));
idx_fp_range1=zeros(1,param.nbr_fp);
idx_fp_range2=zeros(1,param.nbr_fp);
st=1;
for i=1:param.nbr_fp
    vec_idx_fp(st:max_len_fp(i)+st-1)=find(mfp>=fp_vec(i) & mfp<fp_vec(i+1));  
    idx_fp_range1(i)=st;
    idx_fp_range2(i)=max_len_fp(i)+st-1;
    st=st+max_len_fp(i);
end
param.idx_fp_range1=idx_fp_range1;
param.idx_fp_range2=idx_fp_range2;
param.vec_idx_fp=vec_idx_fp;

param.F_dw=fp_vec(1:end-1);
param.F_up=fp_vec(2:end);
param.WF=geomean([fp_vec(2:end);fp_vec(1:end-1)]);
% FPV2=10.^4.547.*param.WF.^(0.938);
FPV2=((1/(10^(-0.665)*1e-6))^(1/0.939))*param.WF.^(1/0.939); %convert back to volume to calculate sinking rate with same dimensions
% param.sink_fp=(10^-2.03.*FPV2.^0.698); %find sinking rates
param.sink_fp=10^(-1.214).*FPV2.^0.513; %changed by the one in "small 1975" paper
% param.sink_fp(param.sink_fp>500)=500;



% param.m_coef=10^-3;%1e-4;5
% 
% m_cte=param.m_coef.*param.Wvec'.^(-1/4);%5e-4;%2e-3;%1e-3;%1e-3;%0.005;

%feeding kernels
param.beta_P=500;%500; %Protists
param.sigma_P=1;
param.beta_act=10000;%1/(10^-3.66); %active feeders
param.sigma_act=1.5;
param.beta_pass=100;%1/10^-2.35; %passive feeders
param.sigma_pass=1; 

%this is the function that defines the thershold for mortality by HTL
% no_small=tanh(param.Wvec');
mmax=max(param.W(:));
beta=param.beta_act;
mshift=max(param.W(:))./beta;
sigma=param.sigma_act;

%--------------------------------
%passive feeders tau
PL=0.532*param.Wvec.^0.32;
SR=1.801*PL-0.695;
sw=10^0.38.*PL.^0.93;
ratio=SR./sw;
tau=ratio;
tau(ratio<0)=0;
tau(param.ind_act)=1;
param.tau=tau;
%this function is to have a reduced mortality for the passive feeders
% p=1/5+tau(param.ind_pass)*(1-1/5); %!!!!!!!!!!!!!!!!!!!
% % p(:)=1;
% 
% phi=exp(-(log(beta*param.W(:)/mmax)).^2./(2*sigma^2));
% phi(param.ind_pass)=phi(param.ind_pass).*p;
% % no_small=phi;
% no_small=phi*0.7+0.3; 
% no_small(param.Wvec>mshift)=1;

p=1/3+tau(param.ind_pass)*(1-1/3); %!!!!!!!!!!!!!!!!!!!
phi=exp(-(log(beta*param.W(:)/mmax)).^2./(2*sigma^2)); %prey  pref function
phi(param.ind_pass)=phi(param.ind_pass)./3; %reduce preference for passives
no_small=phi;%*0.7+0.3; 
no_small(param.Wvec>mshift)=1; %say that any larger than mshift has full preference in mortality htl




% param.md_cte=md_cte;
% param.m_cte=m_cte;
param.p=p;
param.no_small=no_small;

%----------------------------------------
%DETRITUS GRID
WD=zeros(1,length(param.W(:))+param.nbr_P);
param.ind_DP=1:param.nbr_P; %index of dead copepods
param.ind_DC=param.nbr_P+1:length(WD); %index of dead protists
WD(param.ind_DP)=param.V;
WD(param.ind_DC)=param.W(:);

param.nbr_D=nbrD; %number of bins
D_vec=logspace(log10(min(WD))-0.1*min(WD), log10(max(WD)+0.1*max(WD)) ,param.nbr_D+1);
max_len_D=zeros(1,length(D_vec)-1);
for  i=1:param.nbr_D
    max_len_D(i)=length(find(WD>=D_vec(i)  & WD<D_vec(i+1))); 
end

vec_idx_D=zeros(1,length(WD));
idx_D_range1=zeros(1,param.nbr_D);
idx_D_range2=zeros(1,param.nbr_D);
st=1;
for i=1:param.nbr_D
    vec_idx_D(st:max_len_D(i)+st-1)=find(WD>=D_vec(i) & WD<D_vec(i+1));  
    idx_D_range1(i)=st;
    idx_D_range2(i)=max_len_D(i)+st-1;
    st=st+max_len_D(i);
end
param.idx_D_range1=idx_D_range1;
param.idx_D_range2=idx_D_range2;
param.vec_idx_D=vec_idx_D;

param.D_dw=D_vec(1:end-1);
param.D_up=D_vec(2:end);
param.WD=geomean([D_vec(2:end);D_vec(1:end-1)]);

% param.sink_D=20*param.WD.^(0.3); %find sinking rates
param.sink_D=32.*param.WD.^(0.19);

param.m_coef=10^-4;%1e-4;5
% % % % % % 
% % % % % % I_act=1.37*param.Wvec.^(-0.25); %Ingestion rate
% % % % % % k_act=(0.16*param.Wvec.^(-0.25)); %Metabolic cost
% % % % % % clearance_act=(0.011*param.Wvec.^(-0.25)); %Clearance rate
% % % % % % d_act=(param.m_coef.*param.Wvec.^(-0.25)); %mortality rate 0.006
% % % % % % 
% % % % % % I_pass=0.40*param.Wvec.^(-0.25); %Ingestion rate non-calanoid
% % % % % % % I_pass=0.24*param.Wvec.^(-0.25);
% % % % % % k_pass=(0.048*param.Wvec.^(-0.25)); %Metabolic cost thomas non-calanoids 0.0437
% % % % % % clearance_pass=(0.0052*param.Wvec.^(-0.25)); %Clearance rate
% % % % % % d_pass=(param.m_coef.*param.Wvec.^(-0.25));%(0.007.*param.W(param.ind_act).^0); %mortality rate 0.003
% % % % % % 
% % % % % % tau2=tau;
% % % % % % tau2(param.ind_pass)=0;
% % % % % % 
% % % % % % param.I=I_pass+tau2.*(I_act-I_pass);
% % % % % % param.k=k_pass+tau.*(k_act-k_pass);
% % % % % % param.alpha_c=clearance_pass+tau2.*(clearance_act-clearance_pass);
% % % % % % param.d_c=d_pass+tau.*(d_act-d_pass);

I_act=1.37*param.Wvec(param.ind_act).^(-0.25); %Ingestion rate
k_act=(0.16*param.Wvec(param.ind_act).^(-0.25)); %Metabolic cost
clearance_act=(0.011*param.Wvec(param.ind_act).^(-0.25)); %Clearance rate
d_act=(param.m_coef.*param.Wvec(param.ind_act).^(-0.25)); %mortality rate 0.006

I_pass=0.40*param.Wvec(param.ind_pass).^(-0.25); %Ingestion rate non-calanoid
% I_pass=0.24*param.Wvec(param.ind_pass).^(-0.25); %Ingestion rate non-calanoid
% k_pass=(0.048+tau(param.ind_pass).*(0.16-0.048)).*param.Wvec(param.ind_pass).^(-0.25);
k_pass=(0.048*param.Wvec(param.ind_pass).^(-0.25)); %Metabolic cost thomas non-calanoids 0.0437
% k_pass=(0.16*param.Wvec(param.ind_pass).^(-0.25));
clearance_pass=(0.0052*param.Wvec(param.ind_pass).^(-0.25)); %Clearance rate
d_pass=(param.m_coef.*param.Wvec(param.ind_pass).^(-0.25));%(0.007.*param.W(param.ind_act).^0); %mortality rate 0.003

tau2=tau;
tau2(param.ind_pass)=0;

param.I=[I_act; I_pass];
param.k=[k_act; k_pass];
param.alpha_c=[clearance_act; clearance_pass];
param.d_c=[d_act; d_pass];




% Efficiencies
% Assimilation efficiency
% eff_act=0.67*ones(1,length(param.Wvec(param.ind_act)));
% eff_pass=0.67*ones(1,length(param.Wvec(param.ind_pass)));
% param.eff=[eff_act'; eff_pass'];  
param.eff=0.67;

% param.fc=param.k./(param.eff.*param.I); 

% Recruitment efficiencies
% eff_egg_a=0.5;%0.37;
% eff_egg_p=0.5;%0.74;
% 
% % Reproduction efficiencies
% eff_r_a=0.5;%0.42;
% eff_r_p=0.58+tau(param.ind_pass).*(0.42-0.58);
% eff_r_p(:)=0.5;
% 
% %total reproduction efficiency
% rep_eff_act=eff_egg_a*eff_r_a*ones(1,length(param.Wvec(param.ind_act)));
% rep_eff_pass=eff_egg_p.*eff_r_p';%*ones(1,length(param.Wvec(param.ind_pass)));
% param.rep_eff=[rep_eff_act'; rep_eff_pass'];
param.rep_eff=0.25;

% PROTIST --------------------------------------------------------------
cL=21;
AL=0.000914;%0.004864;
invF=linspace(0,1,length(param.V));
% param.alpha_L=((cL*param.V .* AL.*param.V.^(2/3) ./ ( cL*param.V + AL*param.V.^(2/3) )))./param.V; %correct for units of muE to W/m^2
param.alpha_L=(AL*param.V.^(2/3).*(1-exp(-cL.*param.V.^(1/3))))./param.V;
% param.alpha_L=invL.*(AL.*param.V.^(2/3).*(1-exp(-cL2.*param.V.^(1/3))))./param.V;
param.alpha_N=3.75e-5.*param.V.^(-2/3);
param.alpha_F=0.0024.*param.V.^(-1/4);
% param.alpha_F=0.018.*param.V.^(0);

% param.R=0.03*param.V.^(-0.25); %Metabolic cost
param.Qcn=5.6; % C:N ratio

V=2.96e7.*param.V.^1.221;
rho=0.024.*V.^1.1;
Q=0.032.*V.^0.76;
mu_inf=4.7.*V.^-0.26;
param.mu_max=mu_inf.*rho./(mu_inf.*Q + rho);
param.mu_max=param.mu_max + param.mu_max.*0.17;
% param.mu_max(:)=1.5;

param.R=0.2.*param.mu_max; %Metabolic cost


param.m=(10^(-3).*param.V.^-0.25)./(param.V_up./param.V_dw); %mortality rate -4.5 3


% Environment -----------------------------------------------------------
param.kw=0.03;%0.04; %attenuation coefficient of water
param.kchl=1e-5;%3.8e-4; %attenuation coefficient of Chl plankton

% param.Diff_min=0.2;%0.2;%0.13;
% end
param.remin=0.05;
% param.sink=10^-0.03*(1/(25e-9))^0.32*param.WD.^0.32;

% temperature Q10 for physiological rates
param.Q_I=2; %ingestion rate
param.Q_N=1.5; %nitrogen uptake
param.Q_k=2; %respiration
param.Tref=15; %reference temperature

%input of individuals "everything is everywhwre"(?)
param.inputP=(0.001*param.V.^-1).*param.delta_V;
param.inputC=(0.001*param.W(1,:).^-1).*param.deltaC(1,:);
param.inputCa=(0.001*0.1*param.W(end,:).^-1).*param.deltaC(end,:);
param.flow=param.mu_max'.*1e-3; % influx rate
param.flowC=param.I(param.ind_b).*1e-3;%1e-7;
param.flowCa=param.I(param.ind_a).*1e-3;%*1e-3;

%we asume large copepods to migrate, small no
migration=ones(size(param.W(:)));
migration(param.W>1e1)=1;
param.migration=migration;
% param.migration(param.ind_pass)=1;
param.D=10^-2.5;%0.01;

% param.No=100;

deltaC=param.deltaC;
deltaC(end,:)=param.deltaC(param.ind_a-1);
param.deltaC2=deltaC;
deltaratio=param.C_up./param.C_dw;
% deltaratio(end,:)=deltaratio(end-1,:);
param.deltaratio=deltaratio;
%dd_mort_new
% mcoefc=-2;%-3;%3.5;%3;%-4.2;
% param.mmax=(no_small.*10.^(mcoefc).*param.Wvec.^(-0.25))./deltaratio(:);%.*param.C_sp;% 1e-2
% param.mmax=(no_small.*10.^(mcoefc).*param.Wvec.^(-0.25))./deltaC(:);%.*param.C_sp;% 1e-2
param.no_small=no_small;

param.pref_D=zeros(param.nbr_Ctot,param.nbr_D);
param.pref_D(:)=1;

% param.ind_al=param.ind_a(4:5);
% param.Bmin=20;
% param.a_m=1e-3;

v=(param.V/(7.6*1e-7)).^(1/0.819);

% VN=1e-3.*v.^(0.97);
% VN2=(1e-3*(1/(7.6*1e-7)).^(1/0.819)).*m.^((1/0.819)*0.97);
param.VN=2.9617e+04.*param.V.^(1.1844).*5.6*24*1e-6./param.V;
% conv=5.6*24*1e-6./m;

Pc1=10^(-1.45).*v.^(0.3).*24;
Pc2=10^(-0.44).*v.^(-0.2).*24;
param.VL=min(Pc1,Pc2);

param.VF=10^-0.19.*24./1000.*(param.V./1000).^-0.33;

param.mort_coef_protists=0.01;
param.mort_coef_cops=0.01;

end