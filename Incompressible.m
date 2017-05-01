%% 
% 1D Elliptic FD & FV Solver 
% Incompressible Flow Solver - Gaby 24/04/17
%%
% Initialization Parameters, Grid, ...
% Grid Configurations
% |---X---|---X---|...
% Note: # of interfaces = # cells +1
%
clear all
close all
L = 1.0;    % Length of the Reservoir [m]
N = 20;     % Number of Grid Cells

PL = 1; %INPUT BC
PR = 0; %INPUT BC

DX = L/N;   % Grid size
x =     linspace(DX/2, L - DX/2, N);   %Location of Grid centers
xi =    linspace(0, L, N+1); %Location of interfaces

Lambda = zeros(N,1);    %Prealocate Lambda
LmabdaH =zeros(N+1,1); %Prealocate Lambda Harmonic

T = zeros(N+1,1);   %Prealocate Transmissibility

%Lambda(1:N,1) =   1;    %INPUT lambda here  
Lambda(1:N) =   100.*(10.^rand(N,1)); %Rand Het    

LambdaH(1)   = Lambda(1);
LambdaH(N+1) = Lambda(N);

%Compute harmonic averages
%LambdaH(i) = (2*Lambda(i)*Lambda(i-1))/(Lambda(i)+Lambda(i-1));
LambdaH(2:N)  = (2*Lambda(2:N).*Lambda(1:N-1))./(Lambda(2:N)+Lambda(1:N-1));

%Compute Transmissibility values
%T(i) = LambdaH(i)./(DX^2);

T(1)   = LambdaH(1)./(DX^2/2);
T(N+1) = LambdaH(N+1)./(DX^2/2);
T(2:N) = LambdaH(2:N)./(DX^2);

%Plot Results
% Plot permeability and harmonicly averaged fields
% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,'FontSize',14);
plot(x,Lambda,'Parent',axes1,'Marker','.','LineStyle','-','color',[0 0 0],'DisplayName','\lambda');
hold on;
plot(xi,LambdaH,'Parent',axes1,'Marker','*','LineStyle','-', 'color',[1 0 0],'DisplayName','{\lambda}^H');
legend show;
xlabel('x'); 
ylabel(' {\lambda}  and  {\lambda}^H'); 


%%
% Pressure solver
%
A = zeros(N,N);
p = zeros(N,1);
q = zeros(N,1);

%A(2:N-1,2:N-1) = T(2:N-1) + T(3:N);
for i= 1 : N
    if (i>1) %there is a left neighbor
        %T(i) * (p(i)-p(i-1))
        A(i,i) = T(i);
        A(i,i-1) = -T(i);
    end
    if (i<N) %there is a right neighbor
        % T(i) * (p(i)-p(i+1))
        A(i,i) = A(i,i) + T(i+1);
        A(i,i+1) = -T(i+1);
    end 
end

%Insert BC
i = 1; %Left Boundary
% T(i=1) * (P(i=1) -PL)
A(i,i) = A(i,i) + T(i);
q(i) = q(i) + T(i) * PL; 

i = N; %Right Boundary
% T(i+1) * (P(i=N) -PR)
A(i,i) = A(i,i) + T(i+1);
q(i) = q(i) + T(i+1) * PR; 
 p =A\q;
 

 %% Compute the velocity terms
 U = zeros(N+1,1);
 U(1) = 0;
 U(N+1) = 0;
 for i = 2 : N
     U(i) = -LambdaH(i+1) * (p(i)-p(i-1))/DX;
 end
figure
 plot(xi,U);

