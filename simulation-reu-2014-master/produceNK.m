function produceNK(I0)
I0 = 1/600;

global n m numberOfNKCells immuneDensityFine h numberOfImmuneCells
global immuneDensityCoarse immuneDataArray necroticDensityFine

E_nk = I0*(n*m)^2 - numberOfNKCells;
if E_nk <= 0
    newCells = 0;
else
    newCells = randi(round(1+2*E_nk)); % Number of new cells to be placed
end
if newCells > 0
    for k=1:newCells
        done = 0;
        while done==0
            x = randi([1 n*m]);
            y = randi([1 n*m]);
            x1=round(0.5 + h*(x-0.5)); y1=round(0.5 + h*(y-0.5));
            if immuneDensityFine(x,y) == 0
                done = 1;
                immuneDensityFine(x,y) = 1;
                immuneDensityCoarse(x1,y1) = immuneDensityCoarse(x1,y1) + 1;
                immuneDataArray(numberOfImmuneCells+1,:) = [-2 x y x1 y1 0 0 0];
                numberOfImmuneCells=numberOfImmuneCells+1;
                numberOfNKCells=numberOfNKCells+1;
                if necroticDensityFine(x,y) > 0
                    necroticDensityFine(x,y) = necroticDensityFine(x,y) - 1;               
                end
            end
        end
    end
end