function [ nmax ] = Getnmax( ER, UR )
%GETNMAX Method calculates the nmax for a materials array for 1D
%   Detailed explanation goes here

  if(length(ER) ~= length(UR))
    MException('ArraySize', 'Materials arrays are not the same size');
  end

  n = zeros([1 size(ER)]);
  n = sqrt(ER.*UR);  % Calculate refractive index of each material
  nmax = max(n);
end

