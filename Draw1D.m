function h = Draw1D( ER, E, H, dz )
   persistent Color;
   
   %Initialize
   ze=[0:length(E)-1]*dz;
   zh = ze + dz/2;
   
   if isempty(Color)
     % Just inverse our Permitivity to get grayscale value
     ER = ER - min(ER(:));
     
     Color = (ER/max(ER(:)));

     for i = 1 : length(Color)
       if Color(i) < .15 && Color(i) > 0
         Color(i) = 0.15;
       end
     end
     
     Color = abs(Color - 1.10);
     
     for i = 1 : length(Color)
       if Color(i) > 1
         Color(i) = 1;
       end
     end
   end 
   
   % Need to do an initial draw so we can start the hold for plotting.
   cla;   hold on;
   i = 1;
   count = 0;
   prev = 0;
   while i < length(ER)
    i = i + 1;
    
    if(prev == 1)
      prev = ER(i);
      continue;
    end
    
    if(prev == ER(i))
      count = count + 1;
      
    else
      xstart = (i-count)*dz;
      xend = xstart + count*dz;
      
      x = [ xstart xend xend xstart xstart ];
      y = [ -2 -2 2 2 -2 ];
      fill(x,y,[Color(i-1) Color(i-1) Color(i-1)], 'LineStyle', 'none', 'Marker', 'none');
      
      count = 1;
      prev = ER(i);
    end
   end

   %Plot Fields
   h = plot(ze, E, '-b'); 
   
   if(H ~=-1)
     plot(zh, H, '-r'); 
   end;
   
   hold off;
end

