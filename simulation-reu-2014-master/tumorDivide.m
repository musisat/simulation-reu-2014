function tumorDivide(x,y)

global tumorDensityFine tumorDensityCoarse tumorDataArray h numberOfTumorCells n m
global immuneDensityFine hostDensityFine hostDensityCoarse necroticDensityFine tmax

Tmin = tmax;
i_min = 0; % Location with minimum tumor count
j_min = 0;
space = 0;
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
        if tumorDensityFine(x2+i,y2+j)==0 && immuneDensityFine(x2+i,y2+j)==0
            space = space + 1;
        end
    end
end
if space > 0 % If there is a free grid space, choose a random one
    done = 0;
    while done==0
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
        if tumorDensityFine(x2+i,y2+j)==0 && immuneDensityFine(x2+i,y2+j)==0
            done = 1;
            tumorDensityFine(x2+i,y2+j)=1;
            x1_new=round(0.5 + h*(x2+i-0.5)); y1_new=round(0.5 + h*(y2+j-0.5));
            tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
            tumorDataArray(numberOfTumorCells+1,1:5)=[-1 x2+i y2+j x1_new y1_new];
            numberOfTumorCells = numberOfTumorCells + 1;
            if hostDensityFine(x2+i,y2+j) > 0
                hostDensityFine(x2+i,y2+j)=0;
                hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
            elseif necroticDensityFine(x2+i,y2+j) > 0 % necrotic cell
                necroticDensityFine(x2+i,y2+j)=necroticDensityFine(x2+i,y2+j)-1;
            end
        end
    end
 else % Find location with fewest tumor cells
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
                Tcount=tumorDensityFine(x2+i,y2+j);
                if Tcount < Tmin
                    Tmin = Tcount;
                    i_min = i;
                    j_min = j;
                end
            end
        end
        if Tmin<tmax
            if x+i_min<1
                x2=n*m+1;
            elseif x+i_min>n*m
                x2=2;
            else
                x2=x;
            end
            if y+j_min<1
                y2=y+1;
            elseif y+j_min>n*m
                y2=y-1;
            else
                y2=y;
            end
            tumorDensityFine(x2+i_min,y2+j_min)=tumorDensityFine(x2+i_min,y2+j_min)+1;
            x1_new=round(0.5 + h*(x2+i_min-0.5)); y1_new=round(0.5 + h*(y2+j_min-0.5));
            tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
            tumorDataArray(numberOfTumorCells+1,1:5)=[-1 x2+i_min y2+j_min x1_new y1_new];
            numberOfTumorCells = numberOfTumorCells + 1;
            if hostDensityFine(x2+i,y2+j) > 0
                hostDensityFine(x2+i,y2+j)=0;
                hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
            elseif necroticDensityFine(x2+i,y2+j) > 0
                necroticDensityFine(x2+i,y2+j)=necroticDensityFine(x2+i,y2+j)-1;
            end
        end
end



end

