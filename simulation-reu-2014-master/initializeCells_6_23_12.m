function [I0] = initializeCells(radius) 
global n m h tumorDensityCoarse numberOfTumorCells
global immuneDensityFine immuneDensityCoarse numberOfImmuneCells 
global hostDensityFine hostDensityCoarse immuneDataArray tumorDataArray
global killLimit fractionCTL numberOfNKCells tdf highestID
center=round(n*m/2); 
for k=1:numberOfTumorCells
   highestID = highestID + 1;
   done=0; 
   while done==0 
    x=randi([center-radius,center+radius]);
    y=randi([center-radius,center+radius]);
    if (x-center)^2 + (y-center)^2 <= radius^2 
       done=1;
    end
   end
   x1=round(0.5 + h*(x-0.5)); y1=round(0.5 + h*(y-0.5));
   %tumorDensityFine(x,y)=tumorDensityFine(x,y)+1;
   c = tdf(x,y).value;
   tdf(x,y).value = c+1;
   tdf(x,y).id(c+1) = highestID;
   tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)+1;
   tumorDataArray(k,1:5)=[-1,x,y,x1,y1];  % column one = live status (-1 = alive)
   if hostDensityFine(x,y)>0 
       hostDensityFine(x,y)=0; 
       hostDensityCoarse(x1,y1)=hostDensityCoarse(x1,y1)-1;
   end
end

CTLcount = 0;
for k=1:numberOfImmuneCells 
    type=-2; % -2 = NK
    r = rand(); 
    if r < fractionCTL         
        type=-1; % -1 = CTL 
        CTLcount = CTLcount +1;
    end
    done=0; 
    while done==0
       x=randi([2,n*m-1]);   y=randi([2,n*m-1]);
       if tdf(x,y).value==0  
          done=1;    
       end
    end 
    x1=round(0.5 + h*(x-0.5)); y1=round(0.5 + h*(y-0.5));
    immuneDensityFine(x,y)=tdf(x,y).value+1;
    immuneDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)+1;
    immuneDataArray(k,1:5)=[type,x,y,x1,y1]; 
    if type==-1
        immuneDataArray(k,8)=killLimit;
    end
    
end
numberOfNKCells = numberOfImmuneCells - CTLcount;
I0 = numberOfNKCells/(n*m)^2;

end