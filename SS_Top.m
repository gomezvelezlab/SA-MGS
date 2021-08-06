
% Estimation of the Fourier surface and numerical parameters for the
% spectral solutions of groundwater flow systems based on dominant
% frequencies

% This function is an adaptation from Spectop (Wörman et al., 2006)

% Reference: Wörman, A., Packman, A. I., Marklund, L., Harvey, J. W., & Stone, S. H. (2006).
% Exact three-dimensional spectral solution to surface-groundwater interactions with arbitrary 
% surface topography. Geophysical Research Letters, 33(7), 2–5. https://doi.org/10.1029/2006GL025747

% This codes uses a different harmonic function used by Wörman et al., 2006.
% In addition, the Fourier Coefficients are estimated from dominant
% frequencies sampled from "SS_DomFreq". The Regression analysis is
% perfomed using fitlinear which provides more stable parameters avoiding
% ill-conditioned problems in large matrixes.

%   Copyright:  Gabriel Perez
%   Gomez-Velez Lab Group. Vanderbilt University

%   email:   gabriel.perez.mesa@vanderbilt.edu
%	06  Aug	2021,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%	Perez. G., et al., (2021)  Identification of Characteristic Spatial Scales to Improve the Performance 
%   of Analytical Spectral Solutions to the Groundwater Flow Equation (Submitted to WRR)

function[paramhat,kx,ky,h,mean_error,av]=SS_Top(X,Y,Z,Vector_wavelenghts,Ridge_Reg)
    
    % Input arguments

    % X:                                                                    1:n vector containing the x-location in the Z matrix
    % Y:                                                                    1:m vector containing the x-location in the Z matrix
    % Z:                                                                    nxm Matrix containing the DEM in meters
    % Vector_wavelenghts:                                                   Nx2 Matrix with the Wavelenght couples. (This is obtained from the output of SS_Top as [1./fx_dom 1./fy_dom]) 
    % Ridge_Reg:                                                            1: Uses Ridge_Reg to estimate Fourier coeffient (Recommended). 0: Uses Direct method
  
    % Output arguments

    % paramhat:                                                             Amplitude coefficients  with  dimension	(1:Nx2)
    % kx:                                                                   Wave  number,	2*pi/lambda,  in  x-direction with  dimension	(1:N)
    % ky:                                                                   Wave  number,	2*pi/lambda,  in  y-direction with  dimension	(1:N)
    % h:                                                                    Fourier  fitted  surface  that  corresponds  to  the  coordinates X and Y
    % mean_error:                                                           Mean  error  relative  to  the  difference 	(max(Z) 	-  min(Z)).
    % av:                                                                   Average Elevation

% Check Dimension
if  size(X,2)<size(X,1),   X=X';end
if  size(Y,2)<size(Y,1),   Y=Y';end

% Definition  of  lambda-spectrum  included  in  Fourier  series 
av=mean(mean(Z)); %Average  of  Z
N=size(Vector_wavelenghts,1);
% sets of wave  numbers
for  i=1:N
    kx(1,i)=2*pi/Vector_wavelenghts(i,1);%in x-direction
    ky(1,i)=2*pi/Vector_wavelenghts(i,2);%in y-direction
end

% Calculation  of  coefficient  matrix
NX=size(X,2); 
NY=size(Y,2);
if Ridge_Reg==1
        for  j=1:NY
            vect((j-1)*NX+1:(j-1)*NX+NX)=((Z(j,1:NX)-av));
            coef_sin((j-1)*NX+1:(j-1)*NX+NX,1:size(kx,2))=(sin(transpose(X(1:NX))*kx(1:size(kx,2))+Y(j)*ones(NX,1)*ky(1:size(kx,2)))); 
            coef_cos((j-1)*NX+1:(j-1)*NX+NX,1:size(kx,2))=(cos(transpose(X(1:NX))*kx(1:size(kx,2))+Y(j)*ones(NX,1)*ky(1:size(kx,2)))); 
        end
    coef_new=[coef_sin coef_cos]; 
    clear coef_sin coef_cos
    Mdl = fitrlinear(coef_new,vect','Regularization','ridge','solver','lbfgs','learner','svm','FitBias',false);
    clear coef_new vect
    paramhat=Mdl.Beta;
elseif Ridge_Reg==0
        for  j=1:NY
            vect((j-1)*NX+1:(j-1)*NX+NX)=((Z(j,1:NX)-av));
            coef_sin((j-1)*NX+1:(j-1)*NX+NX,1:size(kx,2))=(sin(transpose(X(1:NX))*kx(1:size(kx,2))+Y(j)*ones(NX,1)*ky(1:size(kx,2)))); 
            coef_cos((j-1)*NX+1:(j-1)*NX+NX,1:size(kx,2))=(cos(transpose(X(1:NX))*kx(1:size(kx,2))+Y(j)*ones(NX,1)*ky(1:size(kx,2)))); 
        end
    coef_new=[coef_sin coef_cos]; 
    clear coef_sin coef_cos
    paramhat=pinv(coef_new)*vect'; 
    clear coef_new vect
    %paramhat_new=((coef_new')*coef_new)\((coef_new')*vect'); 

end

% Calculates  surface  and  error  estimates 
for  j=1:NY
    X2D(j,1:NX)=X(1:NX); 
end
for  i=1:NX
    Y2D(1:NY,i)=transpose(Y(1:NY)); 
end

h=av*ones(NY,NX); 
for  p=1:size(paramhat,1)/2 
    h=h+paramhat(p)*sin(kx(p)*X2D+ky(p)*Y2D)+paramhat((size(paramhat,1)/2)+p)*cos(kx(p)*X2D+ky(p)*Y2D);
end 
clear  X2D  Y2D

aaa=h-Z;
mean_error=mean(mean(abs(aaa)))/abs(min(min(Z))-max(max(Z)));
