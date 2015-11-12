function immuneMigrate(x,y,x1,y1,index,iType)

global immuneDataArray immuneDensityFine immuneDensityCoarse h
global tumorDensityFine n m

if iType==-2 % if NK cell, move to location without immune cells
    i = randi([-1 1]); j = randi([-1 1]);
    if x+i<1
        x2=n*m+1;
    elseif x+i>n*m
        x2=2;
    else
        x2=x;
    end
    if y+j<1
        y2=y+1;
    elseif y+j>n*m
        y2=y-1;
    else
        y2=y;
    end
    done = 0;
    guesses = 0;
    while done==0
        if immuneDensityFine(x2+i,y2+j) == 0
            done = 1;
            immuneDensityFine(x2+i,y2+j) = 1;
            immuneDensityFine(x,y)=immuneDensityFine(x,y)-1;
            x1_new=round(0.5 + h*(x2+i-0.5)); y1_new=round(0.5 + h*(y2+j-0.5));
            immuneDensityCoarse(x1_new,y1_new)=immuneDensityCoarse(x1_new,y1_new)+1;
            immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-1;
            immuneDataArray(index,2:5) = [x2+i y2+j x1_new y1_new];
        elseif guesses > 7 % Allow 8 guesses
            done = 1;
        end
        guesses = guesses + 1;
    end
else % if CTL cell, migrate toward highest concentration of tumor cells
    tumorMax = 0; % highest tumor cell count
    i_max = 0; % location of highest tumor cell count
    j_max = 0;
    for i=-1:1
        if x+i<1
            x2=n*m+1;
        elseif x+i>n*m
            x2=2;
        else
            x2=x;
        end
        for j=-1:1
            if y+j<1
                y2=y+1;
            elseif y+j>n*m
                y2=y-1;
            else
                y2=y;
            end
            if tumorDensityFine(x2+i,y2+j) > tumorMax
                tumorMax = tumorDensityFine(x2+i,y2+j);
                i_max = i;
                j_max = j;
            end
        end
    end
    if tumorMax==0
        i_max = randi([-1 1]); j_max = randi([-1 1]);
    end
    if x+i_max<1
        x2=n*m+1;
    elseif x+i_max>n*m
        x2=2;
    else
        x2=x;
    end
    if y+j_max<1
        y2=y+1;
    elseif y+j_max>n*m
        y2=y-1;
    else
        y2=y;
    end
    immuneDensityFine(x2+i_max,y2+j_max)=immuneDensityFine(x2+i_max,y2+j_max)+1;
    immuneDensityFine(x,y)=immuneDensityFine(x,y)-1;
    x1_new=round(0.5 + h*(x2+i_max-0.5)); y1_new=round(0.5 + h*(y2+j_max-0.5));
    immuneDensityCoarse(x1_new,y1_new)=immuneDensityCoarse(x1_new,y1_new)+1;
    immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-1;
    immuneDataArray(index,2:5) = [x2+i_max y2+j_max x1_new y1_new];
end



end

