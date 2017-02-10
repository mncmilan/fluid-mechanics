clear all 
close all
%% %%%%%%%%% FLOW THROUGH A COMPRESSIVE-EXPANSIVE DUCT %%%%%%%%% %%

%% Initial conditions

Sin = 0.5;                  % Inlet section (m)
Sout = 2;                   % Outlet section(m)
theta = 30*pi/180;          % Wall inclination angle (º)
Lin = 0.5;                  % Inlet wall longitude (m)
Lout = 3;                   % Outlet wall longitude (m)


%% Parameters

Dj = 0.2;                   % Separation between points in X axis (m)

%% Previous calculations

% Physical parameters:
Lsl = abs(Sout-Sin)/(round(tan(theta)*10)/10); % Longitud del tramo de expansión (m)

% Numerical parameters:
Di = Dj*(round(tan(theta)*10)/10);  % Separation between points in Y axis (m)
Nj = round((Lin+Lsl+Lout)/Dj + 1);  % Nodes number in X axis
Ni = round(Sout/Di + 1);            % Nodes number in Y axis
beta = Dj^2/Di^2;                   % Beta parameter from equation(8.108) found in White

phi=ones(Ni,Nj);                    % Definition of matrix phi
phiold = zeros(Ni,Nj);              
phiinicial=0;
phi(Ni,1)=phiinicial;               % Initial value of the lower corner

Vinf = 10;
Vold = zeros(Ni,Nj);
V=ones(Ni,Nj);
V(1,1)=Vinf;
V(1,Nj)=Vinf;

%% Boundary conditions

for j=1:Nj
    phi(1,j)=phi(Ni,1)+(Vinf*Sin);       % Upper wall
    phi(Ni,j)=phi(Ni,1);                 % Lower wall
end
for i=2:(Ni-1)
        phi(i,Nj)=phi(i-1,Nj)-(Vinf*(Sin/Sout)*Di);   % Outlet
        phi(i,1)=phi(i-1,1)-(Vinf*Di);                % Inlet
end
phi(Ni,1)=NaN;
for j=1:(Nj-1)  % lower wall
    if(j<round((Lin+Lsl)/Dj))
    phi(Ni,j)= NaN;
    elseif(j==round((Lin+Lsl)/Dj))
        phi(Ni,j)=0;
    end
end

for j=1:(Nj)      % NaNs definition
    for i=Ni-round(round(Lsl*tan(theta)*10)/(10*Dj)):(Ni-1)
        if(j<=round(Lin/Dj))
            phi(i+1,j)=NaN;
            phi(Ni-round(round(Lsl*tan(theta)*10/(Dj*10))),j)=0;
        elseif(j<=round((Lin+Lsl)/Dj))&&(j==(i*round(round(tan(theta)*10)/10))+round((Lin+Lsl)/Dj)-(Ni*round(round(tan(theta)*10)/10)))    
            phi(i,j)=0;           
        elseif(j<round((Lin+Lsl)/Dj))&&(j<(i*round(round(tan(theta)*10)/10))+round((Lin+Lsl)/Dj)-(Ni*round(round(tan(theta)*10)/10)))     
            phi(i,j)=NaN;
        end
    end
end

if((Lsl-(Sout-Sin)*tan(theta))~= 0)
    for j=1:round((Lsl-(Sout-Sin)*tan(theta))/Dj)
    phi(Ni-round(round(Lsl*tan(theta)*10/(Dj*10))),j)=0;
    end
end


 
% DISCRETIZACIÓN ------------------------------------------------------

n=0;
while (abs(nanmean(nanmean(phi-phiold)))>10E-5)
   phiold=phi;
   for j=2:(Nj-1)
      for i=2:(Ni-1)
        if (j>round(Lin/Dj))&&(j>(i*round(round(tan(theta)*10)/10))+round((Lin+Lsl)/Dj)-(Ni*round(round(tan(theta)*10)/10)))
             phi(i,j)=(1/(2*(1+beta)))*(phiold(i,j+1)+phiold(i,j-1)+(beta*phiold(i+1,j))+(beta*phiold(i-1,j)));
        elseif (j<=round(Lin/Dj))&&(i<Ni-round(round(Lsl*tan(theta)*10)/(Dj*10)))
             phi(i,j)=(1/(2*(1+beta)))*(phiold(i,j+1)+phiold(i,j-1)+(beta*phiold(i+1,j))+(beta*phiold(i-1,j)));
        elseif(i<Ni-round(round(Lsl*tan(theta)*10)/(10*Dj)))&&(j<=round(round(Lsl-(Sout-Sin)*tan(theta))/Dj))
             phi(i,j)=(1/(2*(1+beta)))*(phiold(i,j+1)+phiold(i,j-1)+(beta*phiold(i+1,j))+(beta*phiold(i-1,j)));
        end
      end
   end   
   n=n+1;
end

% SIMETRÍA DE LA MATRIZ -----------------------------------------------

phimedia=zeros(Ni,2*Nj);
phitotal=zeros(2*Ni,2*Nj);
for i=1:Ni                              % Simetría izquierda-derecha 
    n=Nj;
    for j=1:Nj
       phimedia(i,j+Nj)=phi(i,j); 
       phimedia(i,j)=phi(i,n); 
       n=n-1;
    end
end

for j=1:2*Nj                            % Simetría arriba-abajo
    for i=1:Ni
        phitotal(i,j) = phimedia(i,j);
    end
    n=Ni;
    for i=(Ni+1):(2*Ni)
        phitotal(i,j) = - phimedia(n,j);
        n=n-1;
    end
end

% PLOT ----------------------------------------------------------------

% Longitud tubería
x=linspace(-Sout,Sout,2*Ni);
y=linspace(0,(Lin+Lsl+Lout)*2,2*Nj);
[X,Y]=meshgrid(y,x);

% Hexágono
xo=[Lout Lout+Lsl Lsl+Lout+(2*Lin) (2*Lin)+Lout+(2*Lsl) Lsl+Lout+(2*Lin) Lout+Lsl Lout];
yo=[0 Lsl*tan(theta) Lsl*tan(theta) 0 -Lsl*tan(theta) -Lsl*tan(theta) 0];

contour(X,Y,phitotal,40)
axis equal
hold on
plot(xo,yo,'Color','black')

