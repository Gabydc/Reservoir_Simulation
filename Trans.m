% This function computes the tansmissibility, taking harmonic averages of
% the permeability 
function [T,LambdaH] = Trans(Lambda, DX, N)

LmabdaH =zeros(N+1,1); %Prealocate Lambda Harmonic

T = zeros(N+1,1);   %Prealocate Transmissibility



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

end