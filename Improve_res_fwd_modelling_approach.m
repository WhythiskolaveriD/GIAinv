function [fwd_model,RMSE] = Improve_res_fwd_modelling_approach(Model,inv_sol,alpha,frad,maskGIA)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% find size of the inv_sol
[r_inv,c_inv] = size(inv_sol);
%find max deg
lmax = r_inv-1;
deg = 1; % assume grid resolution to be 1 degree, change here if need another resolution
% check if the format is CS, if yes then cobnvert it to SC
if r_inv==c_inv
    inv_sol = cs2sc(inv_sol);
end

% define the Gaussian filter
Weights = gaussian(lmax,500);
filter = repmat(Weights,1,(2*(lmax))+1); % filter trends

%check if modle is CS format or SC format
[r_mod,c_mod] = size(Model);
if r_mod == c_mod
    Model_init = Model(1:lmax+1,1:lmax+1);
    Model = cs2sc(Model);
    
else
    mtemp = sc2cs(Model);
    Model_init = mtemp(1:lmax+1,1:lmax+1);
end
lmodel = r_mod - 1;

% create the tranfer function for coverting potential to GIA VLM

tf = (1.1677*(0:lmax) - 0.5233) * 6378137 ;
tf_matrix = repmat(tf',1,(2*(lmax))+1);

tf_full = (1.1677*(0:lmodel) - 0.5233) * 6378137 ;
tf_matrix_full = repmat(tf_full',1,(2*(lmodel))+1);

Weights = gaussian(lmodel,frad);
filter_full = repmat(Weights,1,(2*(lmodel))+1); % filter trends
    
% iterative formward modelling   
% 
% obs_GIA = zeros(size(Model)); Model_GIA_trunc = zeros(size(Model));
% 
% Model_GIA_trunc (1:lmax+1, (lmodel+1 - lmax): (lmodel+1 + lmax)) = cs2sc(Model_init);
% 
% obs_GIA(1:lmax+1, (lmodel+1 - lmax): (lmodel+1 + lmax)) = inv_sol;
% 
% processed_Truth = (Model_GIA_trunc).*filter_full;
% 
% iter = 0 ; RMSE = 1e-12;
%    
%    while RMSE > alpha
%    
%        delta = (obs_GIA.*filter_full) - processed_Truth;
%        %figure; imagesc(delta); mean(delta(:))
%        delta(1:3,:) = 0;
%        updtd_truth_SH = Model + delta;
%        
%        tempTruth = sc2cs(updtd_truth_SH);
%        Model_GIA_trunc (1:lmax+1, (lmodel+1 - lmax): (lmodel+1 + lmax)) = cs2sc(tempTruth(1:lmax+1,1:lmax+1));
%        processed_Truth = (Model_GIA_trunc).*filter_full;
%        
%        RMSE = sqrt(sum(sum(delta.^2))/(121*121));
%        
%        iter = iter + 1;
%    end





      
   Model_grids = gshs(Model.*tf_matrix_full,'none','cell',180,0,0);
   
   processed_Truth = gshs(cs2sc(Model_init).*tf_matrix.*filter,'none','cell',180,0,0);
   
   inv_sol_grids = gshs(inv_sol.*tf_matrix.*filter,'none','cell',180,0,0);
   
   iter = 0 ; RMSE(iter+1) = 1;
   
   updtd_truth_grids = Model_grids;
   
   while RMSE(iter+1) > alpha
   
       delta = (inv_sol_grids - processed_Truth).*maskGIA;
       %figure; imagesc(delta); mean(delta(:))
       
       updtd_truth_grids = updtd_truth_grids + (delta);% - (ones(180,360).*maskGIA*mean(delta(:))));
       
       updtd_truth_SH = gsha(updtd_truth_grids,'mean','block',lmodel);
       
       updtd_truth_SH(1:3,1:3) = sc2cs(cs2sc(Model_init(1:3,1:3)).*tf_matrix(1:3,1:5));
       
       processed_Truth = gshs(cs2sc(updtd_truth_SH(1:lmax+1,1:lmax+1)).*filter,'none','cell',180,0,0);
       
       diff_term = delta;% - (ones(180,360).*maskGIA*mean(delta(:)));
       
       iter = iter + 1;
       
       RMSE(iter + 1) = sqrt(mean(diff_term(:).^2));
       
   end
    
   fwd_model = updtd_truth_SH;
   iter

end

