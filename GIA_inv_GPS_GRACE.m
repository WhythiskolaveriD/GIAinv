function [GIA_svd,Mask_Syn_Data,tf_GIA] = GIA_inv_GPS_GRACE(GRACEtrend,GPStrend,loc_GPS_data,sw,Purc_or_Wahr,prior_model,v_or_e,lmax)
%  GIA_inv_GPS_GRACE implements a novel GIA inversion method that uses contemporary EO data (doi: xxxxx)
%% Input:
% GRACEtrend: Synthetic GRACE data
% GPStrend: nx1 vector containing linear trends in GPS data at n locations. The corodicnates of these locations is provided in the next input, loc_GPS
% loc_GPS: a n x 2 vector with first column containing latitutes and second column containing longitudes
% sw: The grid space between synthetic data sets, we recommend a value of 5
% Purc_or_Wahr: '1' if you want to use relation between VLM and geopotential as provided by Purcell et al 2011, and '2' when you you want to use the relation from Wahr et al 2000 relation.
% prior_model: A GIA model that could be used as a prior and that will help in generating synthetic GNSS dataset.
% v_or_e: '1' if you want to estimate GIA, '2' if you want to estimate PDSMC trends.
% lmax: maximum degree upto which you want to estimate GIA or PDSMC SH coefficients. The maximum value for lmax is dependent of sw, and for sw of 5 lmax should be no more than 40. We recommend 35.

%% Output
% GIA_svd: GIA field obtained without prior
% Mask_Syn_Data: grid map of locations that have been augmented with the synthetic data.
% tf_GIA: transfer function to convert GIA geopotential to VLM.

%% External files/functions used in this function
% gshs, gaussian, constants. These files are available in software bundles available at https://www.gis.uni-stuttgart.de/en/research/downloads/datadrivencorrectionbundle/ (doi:10.1002/2017WR021150)
% or at figshare () doi: xxxxx

%% example %---------------------------------------------------------%
% load('ICE_6GD_stokes.mat')
% load('Caron_stokes.mat')
% load('GRACE_trend_data.mat')
% load('GPS_trend_vec.mat')
% load('GPSloc_vec.mat')

% when using Purcell relation, ICE 6G model as prior (for computing suynthetic data), and solving for GIA
% [GIA_svd_ice6G, ~, tfGIA_purc] = GIA_inv_GPS_GRACE(GRACE_trend,GPStrend,loc_GPS,5,1,CS_ICE_6G_D,1,35);

% same as above but using Caron GIA instead of ICE 6G model as prior
% [GIA_svd_Caron,  ~, tfGIA_purc] = GIA_inv_GPS_GRACE(GRACE_trend,GPStrend,loc_GPS,5,1,Stokes_caron,1,35); when using Caron GIA model as prior (for computing suynthetic data)

% same as above but using Wahr relation instead of Purcell
% [GIA_svd_Caron, ~, tfGIA_wahr] = GIA_inv_GPS_GRACE(GRACE_trend,GPStrend,loc_GPS,5,2,Stokes_caron,1,35); when using ICE 6G model as prior (for computing suynthetic data)

% when using Purcell relation, ICE 6G model as prior (for computing suynthetic data), and solving for PDSMC
% [PDSMC_svd_ice6G, ~, tfGIA_purc] = GIA_inv_GPS_GRACE(GRACE_trend,GPStrend,loc_GPS,5,1,CS_ICE_6G_D,2,35); when using ICE 6G model as prior (for computing suynthetic data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author
% written by Bramha Dutt Vishwakarma, 20th Feb 2021.
% School of Geographical Sciences, University of Bristol.
%% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.

%    see <http://www.gnu.org/licenses/>
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the following paper when using this code:
% Vishwakarma, B.D, Ziegler Y., Bamber, J.L., and Royston S., Separating
% GIA froms urface mass change, 2021, Journal, volume, issue, doi: (to be updated when the article is published)
% main body of the code:

%% check the GIA field and prepare it
[rgm,cgm] = size(prior_model);
GIA_model = prior_model;
if rgm == cgm
    GIACS = cs2sc(GIA_model(1:lmax+1,1:lmax+1));
else
    GIACS = sc2cs(GIA_model);
    GIACS = cs2sc(GIACS(1:lmax+1,1:lmax+1));
end
%% Gaussian filter matrix 

Weights = gaussian(lmax,400);
filter = repmat(Weights,1,(2*(lmax))+1); % filter trends
%% check GRACE data and prepare it
[rGR,cGR] = size(GRACEtrend);
if rGR == cGR
GRACE_tr = GRACEtrend(1:lmax + 1,1:lmax + 1); %.*filter) in m
else
    error('GRACE data not in right format, see input section of function description')
end

%for zero zone constraint
GRACEmGIA = sc2cs(cs2sc(GRACEtrend(1:lmax + 1,1:lmax + 1) - GIACS(1:lmax + 1,1:lmax + 1)));%.*filter);
GRACE_vec =[]; GRACEmGIA_vec =[];
for Gi = 1:lmax+1
    GRACE_vec = [GRACE_vec;GRACE_tr(Gi:lmax+1,Gi)];
    GRACEmGIA_vec = [GRACEmGIA_vec;GRACEmGIA(Gi:lmax+1,Gi)];
    if Gi>1
        GRACE_vec = [GRACE_vec;(GRACE_tr(Gi-1,Gi:lmax+1))'];
        GRACEmGIA_vec = [GRACEmGIA_vec;(GRACEmGIA(Gi-1,Gi:lmax+1))'];
    end
end

%% prepare synthetic GNSS data from GRACE and prior GIA model

GRACEtr_e = gshs(GRACEmGIA,'height','cell',180,0,0);
if Purc_or_Wahr == 1
    tf = (1.1677*(0:lmax) - 0.5233) * 6378137 ;
    filter = repmat(tf',1,(2*(lmax))+1);
    GIA_v = gshs(GIACS.*filter,'none','cell',180,0,0);
elseif Purc_or_Wahr == 2
    tf = (2.*(0:lmax)+1) * 6378137 / 2 ;
    filter = repmat(tf',1,(2*(lmax))+1);
    GIA_v = gshs(GIACS.*filter,'none','cell',180,0,0);
end

GPS_eq = GRACEtr_e/1000 + GIA_v; % in meters

%% transfer functions that relate elastic or viscous VLM with geopotential and EWH with geopotential
constants;
[klGIA,hlGIA,~] = lovenrprem(0:lmax,'CM');
fPurc = (1.1677*(0:lmax) - 0.5233);
if Purc_or_Wahr == 1
    if v_or_e == 1
        tf_GIA =  ( fPurc - (hlGIA ./ (1+klGIA)) ) * ae;% in m * 1000 ; % in mm
        tf_GRACE =  (hlGIA./(1+klGIA)) * ae;
    elseif v_or_e == 2
        tf_GIA =  ( (hlGIA ./ (1+klGIA)) -fPurc) * ae;%
        tf_GRACE =  fPurc * ae;
    else
        error('input not valid for v_or_e');
    end
elseif Purc_or_Wahr == 2
    if v_or_e == 1
        tf_GIA =  ( ((2*(0:lmax)+1).*(1+klGIA) - (2*hlGIA) )./ ((1+klGIA)*2)) * ae;
        tf_GRACE =  (hlGIA./(1+klGIA)) * ae;
    elseif v_or_e == 2
        tf_GIA =  ( ((2*hlGIA) - (2*(0:lmax)+1).*(1+klGIA))./ ((1+klGIA)*2)) * ae;
        tf_GRACE =  ((2*(0:lmax)+1)./2) * ae;
    else
        error('input not valid for v_or_e');
    end
else
    error('input not valid for Purc_or_Wahr');
end

tf_GIA(1:2) = 0; tf_GRACE(2) = 0;

%% use real GNSS data to figure out those 1 degree grid cells that have a GNSS station and then take an average if more than 1 stations are available
ltGPS = 90-loc_GPS_data(:,1); % GPS lat to co-lat
lnGPS = loc_GPS_data(:,2); % GPS longitude
[nl,~] = size(lnGPS);
for li = 1:nl
    if lnGPS(li,1)<0
        lnGPS(li,1)= 361 + lnGPS(li,1);
    end
end

deg = 1;
ln = deg/2:deg:360; lt = deg/2:deg:180;
gp = 1 ;

for ln_s = 1:deg:360
    for lt_s = 1:deg:180
        [clnsearch,~]  = find(lnGPS < ln(ln_s)+(deg/2) & lnGPS > ln(ln_s)-(deg/2) & ltGPS < lt(lt_s)+(deg/2) & ltGPS > lt(lt_s)-(deg/2));
        if ~isempty(clnsearch)
            latGPS(gp) = lt_s-0.5;
            lonGPS(gp) = ln_s-0.5;
            
            GPStr(1,gp) = nanmean(GPStrend(clnsearch));
            
            gp = gp+1;
        end
    end
end
GPS_dl(1:180,1:360) = 0;
for i = 1: gp-1
    GPS_dl(ceil(latGPS(i)),ceil(lonGPS(i))) = 1;
end

thetavec = pi*(latGPS)/180;
lambdavec = pi*(lonGPS)/180;
%% searching for swxsw degree boxes with no data whatsoever
counti =1; counte = 1;
Mask_Data(1:180,1:360) = 0;
Mask_Syn_Data(1:180,1:360) = 0;
for ri = 1:sw:180
    for ci = 1:sw:360
        searchGPS = sum(sum(GPS_dl(ri:ri+sw-1,ci:ci+sw-1)));
        if searchGPS >=1
            Mask_Data (ri:ri+sw-1,ci:ci+sw-1) = 1;
            counti = counti + 1;
        else
            Mask_Data (ri + floor(sw/2), ci + floor(sw/2)) = 0;
            Mask_Syn_Data (ri + floor(sw/2), ci + floor(sw/2)) = 1;
            GPS_vs(counte,1) = mean(mean(GPS_eq(ri:ri+sw-1,ci:ci+sw-1)));
            theta_e(counte,1) = pi*(ri + floor(sw/2) - deg/2)/180;
            lambda_e(counte,1) = pi*(ci + floor(sw/2) - deg/2)/180;
            counte = counte + 1;
        end
    end
end

%% design matrix for GPS and GRACE 
for loc = 1:length(thetavec)
    c_el = 1;
    for ac = 0 : lmax
        if ac == 0
            [Plm] = plm(ac:lmax, 0, thetavec(loc));
            Clambda = cos(0*lambdavec(loc));
            Amat(loc,c_el: c_el + lmax) = Plm*Clambda.*tf_GIA;
            %Amat_GIA(loc,c_el: c_el + lmax) = Plm*Clambda.*tf_GIAfld;
            Amat_GRACE(loc, c_el : c_el + lmax) = Plm*Clambda.*tf_GRACE;
            c_el = lmax + 1;
            
        elseif ac>0
            [Plm] = plm(ac:lmax, ac, thetavec(loc));
            Clambda = cos(ac*lambdavec(loc));
            Slambda = sin(ac*lambdavec(loc));
            Amat(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GIA(ac+1:lmax+1)*Clambda;
            %Amat_GIA(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GIAfld(ac+1:lmax+1)*Clambda;
            Amat_GRACE(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GRACE(ac+1:lmax+1)*Clambda;
            c_el = c_el + length(Plm);
            
            Amat(loc,c_el +1: c_el + length(Plm)) = Plm.*tf_GIA(ac+1:lmax+1)*Slambda;
            %Amat_GIA(loc,c_el +1: c_el + length(Plm)) = Plm.*tf_GIAfld(ac+1:lmax+1)*Slambda;
            Amat_GRACE(loc,c_el +1 : c_el + length(Plm)) = Plm.*tf_GRACE(ac+1:lmax+1)*Slambda;
            c_el = c_el + length(Plm);
        end
    end
end
%% GPS synthetic data to fill empty space and corresponding part of the design matrix
for loc = 1:length(theta_e)
    c_el = 1;
    for ac = 0 : lmax
        if ac == 0
            [Plm] = plm(ac:lmax, 0, theta_e(loc));
            Clambda = cos(0*lambda_e(loc));
            Amat_e(loc,c_el: c_el + lmax) = Plm*Clambda.*tf_GIA;
            %Amat_GIA(loc,c_el: c_el + lmax) = Plm*Clambda.*tf_GIAfld;
            Amat_GRACE_e(loc, c_el : c_el + lmax) = Plm*Clambda.*tf_GRACE;
            c_el = lmax + 1;
            
        elseif ac>0
            [Plm] = plm(ac:lmax, ac, theta_e(loc));
            Clambda = cos(ac*lambda_e(loc));
            Slambda = sin(ac*lambda_e(loc));
            Amat_e(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GIA(ac+1:lmax+1)*Clambda;
            %Amat_GIA(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GIAfld(ac+1:lmax+1)*Clambda;
            Amat_GRACE_e(loc,c_el+1 : c_el + length(Plm)) = Plm.*tf_GRACE(ac+1:lmax+1)*Clambda;
            c_el = c_el + length(Plm);
            
            Amat_e(loc,c_el +1: c_el + length(Plm)) = Plm.*tf_GIA(ac+1:lmax+1)*Slambda;
            %Amat_GIA(loc,c_el +1: c_el + length(Plm)) = Plm.*tf_GIAfld(ac+1:lmax+1)*Slambda;
            Amat_GRACE_e(loc,c_el +1 : c_el + length(Plm)) = Plm.*tf_GRACE(ac+1:lmax+1)*Slambda;
            c_el = c_el + length(Plm);
        end
    end
end

%%

GRACE_eq = Amat_GRACE*GRACE_vec;
GRACE_eq_empty = Amat_GRACE_e*GRACE_vec;


Y = GPStr'/1000 - GRACE_eq;

[~,cA] = size(Amat);%infol = [lmax cA rank(A)];
Acnstrn = zeros(4,cA);
Acnstrn(1,1) = 1; Acnstrn(2,2) = 1; Acnstrn(3,lmax+2) = 1; Acnstrn(4,(2*lmax)+2) = 1;

if v_or_e == 1
Y_empty = [GPS_vs - GRACE_eq_empty;0;0;0;0];
A = [Amat;Amat_e;Acnstrn]; % design matrix
else
    Y_empty = [GPS_vs - GRACE_eq_empty];
    A = [Amat;Amat_e]; % design matrix
end
Y = [Y;Y_empty]; % observation vector

[U,S,V] = svd(A);
iS = pinv(S);
GIAinv = V*iS*U'*Y;% solution or the expectation of the parameters
%% rearrange coefficients into a CS format
GIAfld = GIAinv;
GIAcs(1:lmax+1,1:2*lmax + 1) = 0;
for GVL = 1:lmax+1
    if GVL == 1
        GIAcs(GVL:end,lmax+1) = GIAfld(1:lmax+1); n = lmax+2;
    else
        GIAcs(GVL:end,lmax+GVL) =  GIAfld(n: n + lmax-GVL+1); n = n + lmax + 2 - GVL;
        GIAcs(GVL:end,lmax-GVL+2) = GIAfld(n: n + lmax-GVL+1); n = n + lmax + 2 - GVL;
    end
end
GIA_svd = GIAcs; % output

end

