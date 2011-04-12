function [ output_args ] = Draw2D( xa, ya, ERzz, Ez, NPML, ColorAxis)
%DRAW2D Summary of this function goes here
%   Detailed explanation goes here

  [Nx Ny] = size(Ez);

  cla;   hold on;

  imagesc(xa,ya,Ez);
  caxis([-1*ColorAxis, ColorAxis]);
  
  %Fill in PML  

  if NPML(1)
    x = [xa(1) xa(NPML(1)) xa(NPML(1)) xa(1) xa(1)];
    y = [ya(1) ya(1) ya(Ny) ya(Ny) ya(1)];
    c = 0.5 * [1 1 1];
    fill(x,y,c,'FaceAlpha',0.5);
  end
  if NPML(2)
    x = [xa(Nx-NPML(2)) xa(Nx) xa(Nx) xa(Nx-NPML(2)) xa(Nx-NPML(2))];
    y = [ya(1) ya(1) ya(Ny) ya(Ny) ya(1)];
    c = 0.5 * [1 1 1];
    fill(x,y,c,'FaceAlpha',0.5);
  end
  if NPML(3)
    x = [xa(1) xa(Nx) xa(Nx) xa(1) xa(1)];
    y = [ya(1) ya(1) ya(NPML(3)) ya(NPML(3)) ya(1)];
    c = 0.5 * [1 1 1];
    fill(x,y,c,'FaceAlpha',0.5);
  end
  if NPML(4)
    x = [xa(1) xa(Nx) xa(Nx) xa(1) xa(1)];
    y = [ya(Ny-NPML(4)) ya(Ny-NPML(4)) ya(Ny) ya(Ny) ya(Ny-NPML(4))];
    c = 0.5 * [1 1 1];
    fill(x,y,c,'FaceAlpha',0.5);
  end
  
  
  
  hold off;
  
end

