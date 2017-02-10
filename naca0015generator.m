clear all
close all
%% %%%%%%%%%%%%%%%%%%%%%% Airfoil NACA0015 %%%%%%%%%%%%%%%%%%%%%%% %%

%% Initial parameters
t=0.15; % maximum thickness as a fraction of the chord
c=50; % [mm]chord length
num=100; % number of points

%% Compute airfoil coordinates
X=linspace(0,c,num);
Y=zeros(1,length(X));
i=1;
while i<=length(X)
    Y(i)=(t/0.2)*c*((0.2969*sqrt(X(i)/c))+(-0.1260*(X(i)/c))+((-0.3516*(X(i)/c)^2)+(0.2843*(X(i)/c)^3)+(-0.1015*(X(i)/c)^4)));
    i=i+1;
end
Y(length(X))=0;

% Filling vectors with the coordinates of the lower and upper surface
x=zeros(1,2*length(X)-1);
y=zeros(1,2*length(X)-1);
for i=length(X):length(x)
    x(i)=X(i-length(X)+1);
    y(i)=Y(i-length(X)+1);
end
for i=1:length(X)
    x(i)=X(length(X)-i+1);
    y(i)=-Y(length(X)-i+1);
end
x=x';
y=y';
for i=1:length(x)
    coordinates{i}=[x(i),y(i)];
end
coordinates=cell2mat(coordinates');

%% Plotting airfoil and panels
plot(x,y,'r')
axis equal
title('Airfoil NACA 0015')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')

