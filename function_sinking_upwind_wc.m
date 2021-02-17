function [sinking_stuff, remin_stuff,sinking_flux] = function_sinking_upwind_wc(C,sink_r,deltaz)

% Cb=NaN(nz,nx);
%         for ixz=1:length(C)
%             Cb(nzbox(ixz),nxbox(ixz))=C(ixz);
%         end
% 
Cup=zeros(size(C));
Cup(2:end,:)=C(1:end-1,:);
sinking_stuff = sink_r./deltaz .*(Cup - C);
sinking_flux=sink_r .* C;
% Cbr = sink_r./deltaz .* Cb;

% idx=sub2ind([nz nx],nzbox,nxbox);
% sinking_stuff = Cb2(idx);
remin_stuff = sink_r./deltaz(end) .* C(end);
% Cbcheck=Cb(idx);
% 
% % Check if vector conversion is good
% eq=isequal(Cbcheck,C);
% if eq==0
%    disp('error with vector conversion in sinking upwind function') 
% end


end