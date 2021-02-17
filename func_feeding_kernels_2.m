function [theta_P_P,theta_P_D, theta_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=func_feeding_kernels_2(param)
% create the matrixes for the preference function
% in all: rows are predators, columns are prey
% OUTPUTS:
    % theta_P_P     --> preference of protists for protists
    % theta_P_D     --> preference of protists for detritus
    % theta_P       --> total preference funtion of ptoists (for protists and detritus together)
    % theta_cop_P   --> preference of copepods for protists
    % theta_cop_cop --> preference of copepods for copepods
    % theta_cop_D   --> preference of copepods for detritus
    % theta_cop     --> total preference funtion of copepods (for protists, copepods, and detritus together)

%------------------------------------------------------------------------------------------------------



%calculations of feeding kernels
%rows are predators, columns are prey
theta_P_P=zeros(param.nbr_P,param.nbr_P);

for j=1:param.nbr_P
%arrow function    
theta_P_P(j,:) = sqrt(pi/2)*param.sigma_P* (...
  erf( (log(param.V_dw)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ) ...
  - erf( (log(param.V_up)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ));

end
theta_P_P =theta_P_P./(log(param.V_dw(:))-log(param.V_up(:)));


%feeding preference of generalists on detritus
theta_P_D=zeros(param.nbr_P,param.nbr_D);

for j=1:param.nbr_P
    
theta_P_D(j,:) = sqrt(pi/2)*param.sigma_P* (...
  erf( (log(param.D_dw)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ) ...
  - erf( (log(param.D_up)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ));

end
theta_P_D=theta_P_D./(log(param.D_dw(:)')-log(param.D_up(:)'));
if param.C_sp_pass==0
 theta_P_D(:)=0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

theta_P_D(:)=0;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
theta_P=[theta_P_P, theta_P_D];


%copepods
%it starts with generalists, then species 1 all stages, species 2 all stages etc
%stages are ordered from smallest(top) to adults(bottom)                                                          
pred_vec=param.Wvec;

theta_cop_P=zeros(length(pred_vec),param.nbr_P);

for j=1:length(param.Wvec)
    if param.C_sp_pass>0
        if j<param.ind_pass(1)
            sigma_c=param.sigma_act;
            beta_c=param.beta_act;
        else
%             if param.Wvec(j)<param.oithona_max
%                 sigma_c=param.sigma_pass;
%                 beta_c=param.beta_pass;  
%             else
                sigma_c=param.sigma_pass;
                beta_c=param.beta_pass;  
%             end
        end
    else
        sigma_c=param.sigma_act;
        beta_c=param.beta_act;
    end
    
theta_cop_P(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.V_dw)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.V_up)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));

end
theta_cop_P =theta_cop_P./(log(param.V_dw(:)')-log(param.V_up(:)'));

%prefernce of copepods on copepods
theta_cop_cop=zeros(length(param.W(:)),length(param.W(:)));

for j=1:length(param.W(:))
        if param.C_sp_pass>0
    if j<param.ind_pass(1)
        sigma_c=param.sigma_act;
        beta_c=param.beta_act;
    else
%             if param.Wvec(j)<param.oithona_max
%                 sigma_c=param.sigma_pass;
%                 beta_c=param.beta_pass;  
%             else
                sigma_c=param.sigma_pass;
                beta_c=param.beta_pass;  
%             end             
    end
        else
        sigma_c=param.sigma_act;
        beta_c=param.beta_act; 
            
        end
theta_cop_cop(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.C_dw(:))-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.C_up(:))-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));


end
theta_cop_cop(param.ind_j,param.ind_j) = theta_cop_cop(param.ind_j,param.ind_j)./(log(param.C_dw(param.ind_j))-log(param.C_up(param.ind_j)));
if param.C_sp_pass==0
    theta_cop_cop(param.ind_a,param.ind_j) = theta_cop_cop(param.ind_a,param.ind_j)./(log(param.C_dw(param.ind_j)')-log(param.C_up(param.ind_j)'));
else
    theta_cop_cop(param.ind_a,param.ind_j) = theta_cop_cop(param.ind_a,param.ind_j)./(log(param.C_dw(param.ind_j))-log(param.C_up(param.ind_j)));
end



theta_cop_cop2=theta_cop_cop;

for i=1:length(param.W(:))
            if param.C_sp_pass>0
                if i<param.ind_pass(1)
                    sigma_c=param.sigma_act;
                    beta_c=param.beta_act;
                else
%                     if param.Wvec(j)<param.oithona_max
%                         sigma_c=param.sigma_pass;
%                         beta_c=param.beta_pass;  
%                     else
                        sigma_c=param.sigma_pass;
                        beta_c=param.beta_pass;  
%                     end         
                end
            else
                    sigma_c=param.sigma_act;
                    beta_c=param.beta_act;  
                
            end

theta_cop_cop(i,param.ind_a)=exp(-((log(beta_c*param.W(param.ind_a)./param.W(i))).^2)./(2*sigma_c^2));

end





% lower preference for passive feeding copepods
if param.C_sp_pass>0
% p=tanh(param.W(param.ind_pass)).*2./3+1/3;
% theta_cop_cop(:,param.ind_pass)=p.*theta_cop_cop(:,param.ind_pass); 
% p=tanh(param.W(param.ind_pass)).*2./3+1/3;
% theta_cop_cop(:,param.oithona_idx)=theta_cop_cop(:,param.oithona_idx)./3; 
theta_cop_cop(:,param.ind_pass)=param.p'.*theta_cop_cop(:,param.ind_pass); 
end

theta_cop_F=zeros(length(param.W(:)),param.nbr_fp);

for j=1:length(param.W(:))
    if param.C_sp_pass>0
    if j<param.ind_pass(1)
        sigma_c=param.sigma_act;
        beta_c=param.beta_act;
    else
        sigma_c=param.sigma_pass;
        beta_c=param.beta_pass;              
    end
    else
        sigma_c=param.sigma_pass;
        beta_c=param.beta_pass; 
    end

theta_cop_F(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.F_dw)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.F_up)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));

end
theta_cop_F=theta_cop_F./(log(param.F_dw)-log(param.F_up));
if param.C_sp_pass==0
theta_cop_F(:)=0; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end

% theta_cop_F(:)=0;
theta_cop=[theta_cop_P, theta_cop_cop, theta_cop_F];

% figure
% subplot(2,3,1)
% imagesc(theta_P_P)
% set(gca,'Ydir','normal')
% colorbar
% xlabel('Prey')
% ylabel('Predators')
% title('Preference predation function for Protists')
% 
% subplot(2,3,2)
% imagesc(theta_P_D)
% set(gca,'Ydir','normal')
% colorbar
% xlabel('Prey')
% ylabel('Predators')
% title('Preference predation function for detritus')
% 
% subplot(2,3,4)
% imagesc(theta_cop_P)
% set(gca,'Ydir','normal')
% colorbar
% xlabel('Prey')
% ylabel('Predators')
% title('Preference predation function for copepods on Protists')
% 
% subplot(2,3,5)
% imagesc(theta_cop_cop)
% set(gca,'Ydir','normal')
% colorbar
% xlabel('Prey')
% ylabel('Predators')
% title('Preference predation function for copepods on Copepods')
% 
% subplot(2,3,6)
% imagesc(theta_cop_D)
% set(gca,'Ydir','normal')
% colorbar
% xlabel('Prey')
% ylabel('Predators')
% title('Preference predation function for copepods on detritus')

end