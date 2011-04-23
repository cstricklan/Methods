function [ nmax ] = Getnmax2D( ER, UR )
%GETNMAX2D Method calculates the nmax for a materials array for 2D
%   Detailed explanation goes here


  [eM eN] = size(ER);
  [uM uN] = size(UR);
  
  if(eM ~= uM && eN ~= uN)
    MException('ArraySize', 'Materials arrays are not the same size');
  end

  n = zeros([eM eN]);
  n = sqrt(ER.*UR);  % Calculate refractive index of each material
  nmax = max(max(n));
end




