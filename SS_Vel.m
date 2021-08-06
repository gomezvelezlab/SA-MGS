% Exact Solution to groundwater flow in 3D based on dominant
% frequencies

% This function is an adaptation from Specvel (Wörman et al., 2006)

% Reference: Wörman, A., Packman, A. I., Marklund, L., Harvey, J. W., & Stone, S. H. (2006).
% Exact three-dimensional spectral solution to surface-groundwater interactions with arbitrary 
% surface topography. Geophysical Research Letters, 33(7), 2–5. https://doi.org/10.1029/2006GL025747

% This codes uses a different harmonic function used by Wörman et al., 2006.
% The inputs are provided by SS_top routine.

%   Copyright:  Gabriel Perez
%   Gomez-Velez Lab Group. Vanderbilt University

%   email:   gabriel.perez.mesa@vanderbilt.edu
%	06  Aug	2021,   MATLAB	2019b  version
%	IF  YOU	PUBLISH  WORK  BENEFITING  FROM  THIS  M-FILE,   PLEASE  CITE  IT AS:
%	Perez. G., et al., (2021)  Identification of Characteristic Spatial Scales to Improve the Performance 
%   of Analytical Spectral Solutions to the Groundwater Flow Equation (Submitted to WRR)

function[U,V,W,x,y,z,h]=SS_Vel(X,Y,fact,N_Z,av,kx,ky,paramhat,Cond,poro,z_min,z_max,dp)

% Input arguments

    % X:                                                                    1:n vector containing the x-location in the Z matrix
    % Y:                                                                    1:m vector containing the x-location in the Z matrix
    % fact:                                                                 Multiples of grid step used in X and Y. For computational efficiency,
                                                                          % the original discretisation DX_orig = X(2) - X(1) i s increased as DX_new = fact*DX_orig
    % N_Z:                                                                  Number of grid steps with depth 
    % av:                                                                   Average value of the estimated surface.
    % kx:                                                                   Wave  number,	2*pi/lambda,  in  x-direction with  dimension	(1:N)
    % ky:                                                                   Wave  number,	2*pi/lambda,  in  y-direction with  dimension	(1:N)
    % paramhat:                                                             Amplitude coefficients  with  dimension	(1:Nx2)
    % Cond:                                                                 Hydraulic conductivity (m/s)
    % poro:                                                                 Porosity
    % z_min:                                                                Min Depth in meters in the output
    % z_max:                                                                Max Depth in meters in the output (This can be equal to dp)
    % dp:                                                                   Depth flow domain in meters. Recommended to be at least 2 or 3 times the max domain lenght
  
    % Output arguments

    % U:                                                                   Velocity vector in x -direction with dimensions (1:m,1:n,1:N_Z+1)
    % V:                                                                   Velocity vector in y -direction with dimensions (1:m,1:n,1:N_Z+1)
    % W:                                                                   Velocity vector in z -direction with dimensions (1:m,1:n,1:N_Z+1)
    % x:                                                                   Coordinate in x -direction with dimensions (1:m,1:n,1:N_Z+1)
    % y:                                                                   Coordinate in y -direction with dimensions (1:m,1:n,1:N_Z+1)
    % z:                                                                   Coordinate in z -direction with dimensions (1:m,1:n,1:N_Z+1)
    % h:                                                                   Fourier fitted hydraulic potential that corresponds to the grid (x,y,z) with dimensions (1:m,1:n,1:N_Z+1)

if size(X,2)<size(X,1), X=X';end 
if size(Y,2)<size(Y,1), Y=Y';end
    x_min=min(min(X)); % Coordinates of volume
    x_max=max(max(X));
    y_min=min(min(Y));
    y_max=max(max(Y));
if dp>40000, dp=40000; end
Cond_ef=Cond/poro;

DX=fact*(x_max-x_min)/(size(X,2)-1);
DY=fact*(y_max-y_min)/(size(Y,2)-1);
N_X=floor(size(X,2)/fact)-1;
N_Y=floor(size(Y,2)/fact)-1;
x_max=x_min+N_X*DX;
y_max=y_min+N_Y*DY;

N=size(paramhat,2)/2; % GP
if size(paramhat,1)>size(paramhat,2), N=size(paramhat,1)/2; end 

% Note that in Matlabs plot routines, x corresponds to the columns of
% func(j,i) and y corresponds to the rows
for i=1:N_X+1
    for j=1:N_Y+1
        for k=1:N_Z+1
            x3D(j,i,k)=x_min+(i-1)*(x_max-x_min)/N_X;
            y3D(j,i,k)=y_min+(j-1)*(y_max-y_min)/N_Y;
            z3D(j,i,k)=z_min+(k-1)*(z_max-z_min)/N_Z;
        end
    end
end
x(1:N_X+1)=x3D(1,1:N_X+1,1);
y(1:N_Y+1)=y3D(1:N_Y+1,1,1);
z(1:N_Z+1)=z3D(1,1,1:N_Z+1);
h=av*ones(N_Y+1,N_X+1,N_Z+1);
U=zeros(N_Y+1,N_X+1,N_Z+1);
V=zeros(N_Y+1,N_X+1,N_Z+1);
W=zeros(N_Y+1,N_X+1,N_Z+1);
for p=1:N 
    Umx=Cond_ef*kx(p);
    Umy=Cond_ef*ky(p);
    Umz=Cond_ef*sqrt(kx(p)^2+ky(p)^2);
    nom=exp(z3D*sqrt(kx(p)^2+ky(p)^2))+exp(sqrt(kx(p)^2+ky(p)^2)*(-2*dp-z3D));
    denom=1+exp(-sqrt(kx(p)^2+ky(p)^2)*2*dp);
    F1=nom/denom;
    h=h+F1.*(paramhat(p).*sin(x3D*kx(p)+y3D*ky(p))+paramhat(N+p).*cos(x3D*kx(p)+y3D*ky(p)));
    U=U-Umx*F1.*(paramhat(p).*cos(x3D*kx(p)+y3D*ky(p))-paramhat(N+p).*sin(x3D*kx(p)+y3D*ky(p)));
    V=V-Umy*F1.*(paramhat(p).*cos(x3D*kx(p)+y3D*ky(p))-paramhat(N+p).*sin(x3D*kx(p)+y3D*ky(p)));
    nom=exp(sqrt(kx(p)^2+ky(p)^2)*z3D)-exp(sqrt(kx(p)^2+ky(p)^2)*(-2*dp-z3D));
    denom=1+exp(-sqrt(kx(p)^2+ky(p)^2)*2*dp);
    F2=nom/denom;
    W=W-Umz*F2.*(paramhat(p).*sin(x3D*kx(p)+y3D*ky(p))+paramhat(N+p).*cos(x3D*kx(p)+y3D*ky(p)));
end

