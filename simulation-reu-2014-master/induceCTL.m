function induceCTL(x,y)
% Check the surrounding cells for CTL induction.

global tumorDensityFine immuneDensityFine immuneDensityCoarse thetaL
global immuneDataArray killLimit numberOfImmuneCells h
 
T_near = 0; % Tumor cells in neighborhood
for i=-1:1
    for j=-1:1
        T_near=T_near+tumorDensityFine(x+i,y+j);
    end
end
for i=-1:1
    for j=-1:1
        if tumorDensityFine(x+i,y+j)==0 && immuneDensityFine(x+i,y+j)==0 && ~(i==0 && j==0)
            r = rand();
            P_L=exp(-(thetaL/T_near)^2);
            if r < P_L % Recruit new CTL cell to that location
                x1=round(0.5 + h*(x+i-0.5)); y1=round(0.5 + h*(y+j-0.5));
                immuneDensityFine(x+i,y+j)=1;
                immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1);
                immuneDataArray(numberOfImmuneCells+1,:)=[-1 x+i y+j x1 y1 0 0 killLimit];
                numberOfImmuneCells = numberOfImmuneCells+1;
            end
        end 
    end
end

end

