function [ N_lambda ] = GetNlambda( ER, UR )
%GETNLAMBDA Summary of this function goes here
%   Detailed explanation goes here

x = size(ER);

if (x(1) == 1)
  nmax = Getnmax(ER, UR);
else
  nmax = Getnmax2D(ER, UR);
end

if(nmax < 10)
  N_lambda = 20;
end

if((nmax > 10) && (nmax<=40))
  N_lambda = 30;
end

if((nmax > 40) && (nmax <= 60))
  N_lambda = 60;
end


if(nmax > 60)
  N_lambda = 200;
end


end

