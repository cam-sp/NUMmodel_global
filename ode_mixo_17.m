function [dydt,s_flux_D_diag, s_flux_F_diag, NPP, pred_P_diag, pred_C_on_P_diag, pred_C_on_C_diag, pred_C_on_D_diag, pred_C_on_F_diag,dd_mort_C_diag,com_resp,gamma_a,pred_C_Dtot] =...
    ode_mixo_17(t, y, param,theta_P, theta_cop_P, theta_cop_cop, theta_cop_F, Ir, nbrP,...
    I,alpha_N,mu_max,alpha_F,alpha_c,k,layer,z,deltaz,remin,VN,VF,VL,diags,theta,idxBmort,lat,lon)


if any(isnan(y))
      warning('there are NaNs inside the ODE!!=');
      disp(lat)
      disp(lon)
end  

% to avoid too small time-stepping of the solver
if any(y<1e-50)
   y(y<1e-50)=1e-50; 
end


alpha_c=alpha_c';

N = y(1:layer);
P1 = y(layer+1:layer+layer*param.nbr_P);
C1 = y(layer+layer*param.nbr_P+1:layer+layer*param.nbr_P+layer*param.nbr_Ctot);
F1 = y(layer+layer*param.nbr_P+layer*param.nbr_Ctot+1:end-layer*param.nbr_D);
D1 = y(end-layer*param.nbr_D+1:end);

P=reshape(P1,[],nbrP);
C=reshape(C1,[],param.nbr_Ctot);
F=reshape(F1,[],param.nbr_fp);
D=reshape(D1,[],param.nbr_D);

% if any(C<1e-20)
%     C(C<1e-20)=0;
% end  


% I0=Ir;
% Il=zeros(layer,1);
% Pgrid=cumsum(sum(P,2).*z(1:layer),1);
% for i=1:layer
%     ktot = param.kw+ param.kchl .* Pgrid(i); % attenuation coeficient
%     Il(i)=I0./(ktot.*deltaz(i)).*(1-exp(-ktot.*deltaz(i))); %avg light in the mixed layer
%     Izh=Ir.*exp(-ktot.*z(i)); %light in the depth horizon
%     I0=Izh;
% end
L=Ir;


    E_P_tot=P*theta_P';
    E_C_tot=(theta_cop_P*P')' + (theta_cop_cop*C')' + (param.pref_D*D')' + (theta_cop_F*F')';
% E_C_tot(:)=0; %!!!!!!!!!!!!!!!!!!!!!!!!!

%     J_F=(mu_max.*alpha_F.*E_P_tot)./(alpha_F.*E_P_tot+mu_max);
%     J_N=(mu_max.*param.Qcn.*alpha_N.*N)./(param.Qcn.*alpha_N.*N+mu_max);
%     J_L=(mu_max.*param.alpha_L.*L)./(param.alpha_L.*L+mu_max);
    
    J_F=(VF.*alpha_F.*E_P_tot)./(alpha_F.*E_P_tot+VF);
    J_N=(VN.*param.Qcn.*alpha_N.*N)./(param.Qcn.*alpha_N.*N+VN);
    J_L=(VL.*param.alpha_L.*L)./(param.alpha_L.*L+VL);

    J_R=0.17.*mu_max;

%     mu=zeros(size(J_L));
%     realtheta=zeros(size(J_L));
% %     theta=1;

        muJCtot=J_L+J_F- J_R;
        muJNtot=J_N+ J_F;
        mu=min(muJCtot,muJNtot);
   
    leaks=max(0,J_N - J_L + J_R);
    
    pred_P=(J_F.*P./E_P_tot)*theta_P;
%     pred_P=((1-realtheta).*J_F.*P./E_P_tot)*theta_P;
    
    
%     pred_P_on_D=((alpha_F.*(1-fg).*P)*theta_D); 
   
    %Feeding level of copepods
     F_lvl=(alpha_c'.*E_C_tot)./(alpha_c'.*E_C_tot+I);
% k=param.eff.*I.*(0.1+0.2.*F_lvl);
    %Net energy gain for all stages
    fc=k./(param.eff'.*I);
    nu=param.eff'.*I.*(F_lvl-fc); 

    
    pred_C_on_P=(I.*F_lvl.*C./E_C_tot)*theta_cop_P;
    
    %background mortality of protists
    mortP=mu_max.*param.mort_coef_protists.*P./param.ratio_V; %0.03
    
    % Protists ODE
    dPdt = mu.*P -  mortP.*P - pred_P.*P - pred_C_on_P.*P + param.flow'.*(param.inputP - P); %last term is to ensure a small bckground concentration

 
mmax=(param.no_small'.*I.*param.mort_coef_cops)./param.deltaratio(:)';% 0.015 0.003
dd_mort_C=zeros(size(C));
% dd_mort_C=mmax.*C;
for i=1:length(param.Wvec)
   Bcop=sum(C(:,idxBmort{i}),2);
   dd_mort_C(:,i)=mmax(:,i).*C(:,i).^(0.2).*Bcop.^(1-0.2);    
end

dCdt= - dd_mort_C.* C;

    nu_pos=max(0,nu); %positive growth
    nu_neg=-min(0,nu); %starvation mortalityt
    nu_ca_neg=min(0,nu); %starvation of the adults
    
    pred_C_on_C=(I.*F_lvl.*C./E_C_tot)*theta_cop_cop;
   
%     pred_C_on_C(1,:)=(I(1,:).*F_lvl(1,:).*C(1,:)./E_C_tot(1,:))*theta_cop_cop(:,param.ind_a)
    
    %Predation from copepods on fecal pellets
    pred_C_on_F=(I.*F_lvl.*C./E_C_tot)*theta_cop_F;


    %Predation from copepods on detritus
    pred_C_on_D=(I.*F_lvl.*C./E_C_tot)*param.pref_D;


    %maturation rates
    d_c=pred_C_on_C+nu_neg+dd_mort_C;
    gamma=zeros(size(C)); %for convinience in later calculations we make it in this form
    gamma(:,param.ind_j)=(nu_pos(:,param.ind_j)-d_c(:,param.ind_j))./(1-param.z(param.ind_j)'.^(1-(d_c(:,param.ind_j)./nu_pos(:,param.ind_j))));
    
    
    %ODE copepod
    % Size at birth
    dCdt(:,param.ind_b)=dCdt(:,param.ind_b) + param.rep_eff.* nu_pos(:,param.ind_a).*C(:,param.ind_a)...
                     + nu(:,param.ind_b).*C(:,param.ind_b) ...
                     - gamma(:,param.ind_b).*C(:,param.ind_b) - pred_C_on_C(:,param.ind_b).*C(:,param.ind_b) ...
                     + param.flowC'.*(param.inputC - C(:,param.ind_b));
    
    % In-between size-classes
    dCdt(:,param.ind_rest)=dCdt(:,param.ind_rest) +  gamma(:,param.ind_rest-1).*C(:,param.ind_rest-1) + nu(:,param.ind_rest).*C(:,param.ind_rest) ...
                          - gamma(:,param.ind_rest).*C(:,param.ind_rest)...
                          - pred_C_on_C(:,param.ind_rest).*C(:,param.ind_rest);
                                                
    % Adults                  
    dCdt(:,param.ind_a) = dCdt(:,param.ind_a) + gamma(:,param.ind_a-1).*C(:,param.ind_a-1) + nu_ca_neg(:,param.ind_a).*C(:,param.ind_a) ...
                        - pred_C_on_C(:,param.ind_a).*C(:,param.ind_a);

                    
% Fecal pellet production                
fpp=(1-param.eff').*I.*F_lvl;
FPP=fpp.*C;

dFdt=zeros(size(F)); %vector to fill
for i=1:param.nbr_fp
    range_fp=param.vec_idx_fp(param.idx_fp_range1(i):param.idx_fp_range2(i));
    dFdt(:,i)=sum(FPP(:,range_fp),2);
end
%fecal pellets ode
dFdt = dFdt - remin.*F  - pred_C_on_F.*F;


% Phyto-detritus production                
DP=mortP.*P;%param.m.*P.^2;
DC=(0.1.*dd_mort_C).*C; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DPP=[DP, DC];
dDdt=zeros(size(D)); %vector to fill
for i=1:param.nbr_D
    range_D=param.vec_idx_D(param.idx_D_range1(i):param.idx_D_range2(i));
    dDdt(:,i)=sum(DPP(:,range_D),2);
end
dDdt= dDdt - remin.*D - pred_C_on_D .*D;
  
dNdt=( - sum(J_N .* P, 2) + sum((leaks.*P),2) + remin.*(sum(D,2) + sum(F,2)) + sum(k.*C,2)...
         + sum((1-param.rep_eff).* nu_pos(:,param.ind_a).*C(:,param.ind_a),2) ...
         + sum(0.9.*dd_mort_C.*C,2))./param.Qcn; % in copepods respiration 

                 
            
                if any(isnan(F))
                  disp('there are NaNs before sinking F!!=');

                end  
                
                if any(isnan(D))
                  disp('there are NaNs before sinking D!!=');

                end  


%do stuff sinking   

%----------------------------------
s_flux_D_diag=0;
for i=1:param.nbr_D
    [Dsink, Dremin, s_flux_D] = function_sinking_upwind_wc(D(:,i),param.sink_D(i),deltaz(1:layer));
    dDdt(:,i)=dDdt(:,i) + Dsink;
    dNdt(end)=dNdt(end) + Dremin./param.Qcn;% .*idx_land_interface;
    s_flux_D_diag=s_flux_D_diag+s_flux_D;
end

s_flux_F_diag=0;
for i=1:param.nbr_fp
    [Fsink, Fremin, s_flux_F] = function_sinking_upwind_wc(F(:,i),param.sink_fp(i),deltaz(1:layer));
    dFdt(:,i)=dFdt(:,i) + Fsink;
    dNdt(end)=dNdt(end) + Fremin./param.Qcn;% .*idx_land_interface;
    s_flux_F_diag=s_flux_F_diag+s_flux_F;
end

%----------------------------------    
            
                if any(isnan(F))
                  disp('there are NaNs after sinking F!!=');
                  
                end  
                
                if any(isnan(D))
                  disp('there are NaNs FTER sinking D!!=');

                end            
            
            
dPdt2=reshape(dPdt,[],1);
dCdt2=reshape(dCdt,[],1);
dFdt2=reshape(dFdt,[],1);
dDdt2=reshape(dDdt,[],1);

%check mass balance/ note that a correction should be done for the sinking
%rates, so if we have dinking rates different of 0 we should correct for
%the width of the layer to properly calculate mass balance. Here we can
%check mass balance 
% if t<0.01
%  
% dNtotdt=sum((dNdt + (sum(dPdt,2)+ sum(dCdt,2)+ sum(dFdt,2) +sum(dDdt,2))./param.Qcn).*deltaz(1:layer));
% dNtotdt=dNtotdt + (- sum(sum((param.flowC'.*(param.inputC - C(:,param.ind_b))),2).*deltaz(1:layer))...
%     - sum(sum((param.flow'.*(param.inputP - P)),2).*deltaz(1:layer)))./param.Qcn; %last term is to ensure a small bckground concentration

%     if dNtotdt>1e-10
%      disp('no mass balance')
%     end
% end
    
if diags==0

dydt = [dNdt; dPdt2; dCdt2; dFdt2; dDdt2];
% F_lvl=0; 
% d_c=0;
% nu=0;
gamma_a=0;
NPP=0;
elseif diags==1
dydt =0;
% sinking_det=Dsink
% sinking_fp= 
NPP=min(J_L,J_N).*P;
% GPP=min(realtheta.*J_L,J_N).*P;

pred_P_diag=pred_P.*P;
pred_C_on_P_diag=pred_C_on_P.*P;
pred_C_on_C_diag=pred_C_on_C.*C;
pred_C_on_D_diag=pred_C_on_D.*D; 
pred_C_on_F_diag=pred_C_on_F.*F;
pred_C_Dtot=sum(pred_C_on_D_diag,2)+sum(pred_C_on_F_diag,2);
dd_mort_C_diag=dd_mort_C.*C;

% gamma_a=gamma(:,param.ind_a-1).*C(:,param.ind_a-1).*deltaz(1:layer);
gamma_a=(sum((0.9.*dd_mort_C).*C,2)).*deltaz(1:layer);

com_resp=(sum((1-param.rep_eff).* nu_pos(:,param.ind_a).*C(:,param.ind_a),2) + sum(k.*C,2) + sum(J_R.*P,2)).*deltaz(1:layer);


end

end