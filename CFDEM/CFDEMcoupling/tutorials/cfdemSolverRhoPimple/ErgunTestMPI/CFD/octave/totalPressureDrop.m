%
close all;
clear;
clc;

%====================================%
% simulation data 1
%====================================%
rhoG = 1.18;
path = '../postProcessing/probes/0/p';
data = load(path);
dp_sim = (data(:,2)-data(:,3));
t_sim = data(:,1);

%====================================%
% analytical calculation
%====================================%

%===================
% Ergun Equation
%===================

dp = 0.001;
phip = 1;
epsilon = 0.451335;
Ustart = 0.01;
Uend = 1.0;
deltaU= (Uend-Ustart)/(length(t_sim)-1);
U = Ustart:deltaU:Uend;
L = 0.0156;
muG = 2e-05;

dpErgun= L * (
                150*((1-epsilon)^2/epsilon^3)*((muG.*U)/(phip*dp)^2) 
              +1.75*((1-epsilon)/epsilon^3)*((rhoG.*U.^2)/(phip*dp))
        );

%==================================
% min fluidization velocity in m/s
%==================================

rhoP = 2000;
g = 9.81;

a = 1.75*((1-epsilon)/epsilon^3)*(rhoG/(phip*dp));
b = 150*((1-epsilon)^2/epsilon^3)*(muG/(phip*dp)^2); 
c = -(rhoP-rhoG)*g*(1-epsilon);
sqrtD = (b^2-4.0*a*c)^0.5;
Umf = (-b + sqrtD)/(2.0*a);

ReMF = Umf*dp*rhoG/muG;

dpUmf= L * (
                150*((1-epsilon)^2/epsilon^3)*((muG.*Umf)/(phip*dp)^2) 
              +1.75*((1-epsilon)/epsilon^3)*((rhoG.*Umf.^2)/(phip*dp))
        );

%====================================%
% plot data
%====================================%

fig=figure();
title("Ergun pressure drop vs. simulation")
plot(U,dpErgun,U,dp_sim,[Umf,Uend],dpUmf*ones(1,2))
a=strcat("analytical (Ergun), Umf=",num2str(Umf),", dpUmf=",num2str(dpUmf));
legend(a,"simulation","analyt. deltaP at Umf","location","northwest")
xlabel("velocity in [m/s]")
ylabel("pressure drop [Pa]")
axis([0,Uend,0,dpErgun(length(dpErgun))])
print(fig,"cfdemSolverPiso_ErgunTestMPI.png")
