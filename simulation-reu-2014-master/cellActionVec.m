function cellActionVec(cellType)

global immuneDataArray tdf thetaLD killLimit tumorDataArray thetaL
global immuneDensityFine numberOfNKCells theta_nec tmax tumorDensityCoarse
global numberOfTumorCells theta_div theta_mig n m Mmin Nmin hostDensityCoarse
global hostDensityFine necroticDensityFine h numberOfImmuneCells immuneDensityCoarse

if cellType == 1 % Tumor cells
    action = randi([1 3], [numberOfTumorCells 1]);
    idMig = find(action==1);
    idNec = find(action==2);
    idDiv = find(action==3);

    % Migrate:
    if numel(idMig) > 0
        r1 = rand(size(idMig));
        z1 = tumorDataArray(idMig,6); % Find M nutrient level
        P_mig = 1-exp(-numberOfTumorCells*(z1/theta_mig).^2);
        go1 = r1 < P_mig & tumorDataArray(idMig,1) == -1;
        idTemp1 = find(go1);
        idMigFinal = idMig(idTemp1);
        for k=1:numel(idMigFinal)
            index=idMigFinal(k);
            x=tumorDataArray(index,2);
            y=tumorDataArray(index,3);
            x1=tumorDataArray(index,4);
            y1=tumorDataArray(index,5);
            if x~=0 && y~=0
                Tmin = tmax; % Minimum tumor cell count
                i_min = 0; % Location with minimum tumor count
                j_min = 0;
                count = 1;
                space = zeros(8,2);
                if tdf(x,y).value > 0
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
                            if ~(i==0 && j==0) && tdf(x2+i,y2+j).value==0 && immuneDensityFine(x2+i,y2+j)==0
                                space(count,:) = [x2+i y2+j];
                                count = count + 1;
                            end
                        end
                    end
                    if numel(find(space)) > 0 % If there is a free grid space, choose a random one
                        spaceFinal = space(find(space(:,1)),:);
                        s = size(spaceFinal);
                        ri = randi([1 s(1)]);
                        chosenX = spaceFinal(ri,1);
                        chosenY = spaceFinal(ri,2);
                        tdf(chosenX,chosenY).value = 1;
                        c = tdf(x,y).value;
                        tdf(x,y).id(c) = 0;              
                        tdf(x,y).value = c-1;
                        tdf(chosenX,chosenY).id(1) = index;
                        tdf(chosenX,chosenY).value = 1;
                        x1_new=round(0.5 + h*(chosenX-0.5));
                        y1_new=round(0.5 + h*(chosenY-0.5));
                        tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
                        tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)-(tumorDensityCoarse(x1,y1)>0);
                        tumorDataArray(index,2:5) = [chosenX chosenY x1_new y1_new];
                        if hostDensityFine(chosenX,chosenY) > 0
                            hostDensityFine(chosenX,chosenY)=0;
                            hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
                        elseif necroticDensityFine(chosenX,chosenY) > 0 % necrotic cell
                            necroticDensityFine(chosenX,chosenY)=necroticDensityFine(chosenX,chosenY)-1;
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
                                Tcount=tdf(x2+i,y2+j).value;
                                if Tcount < Tmin
                                    Tmin = Tcount;
                                    i_min = i;
                                    j_min = j;
                                end
                            end
                        end
                        if Tmin < tmax && Tmin > 0
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
                            c = tdf(x2+i_min,y2+j_min).value;
                            tdf(x2+i_min,y2+j_min).value = c+1;
                            tdf(x2+i_min,y2+j_min).id(c+1) = index;
                            d = tdf(x,y).value;
                            tdf(x,y).id(d) = 0;
                            tdf(x,y).value = d-1;
                            x1_new=round(0.5 + h*(x2+i_min-0.5));
                            y1_new=round(0.5 + h*(y2+j_min-0.5));
                            tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
                            tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)-(tumorDensityCoarse(x1,y1)>0);
                            tumorDataArray(index,2:5)=[x2+i_min y2+j_min x1_new y1_new];
                            if hostDensityFine(x2+i,y2+j) > 0
                                hostDensityFine(x2+i,y2+j)=0;
                                hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
                            elseif necroticDensityFine(x2+i,y2+j) > 0
                                necroticDensityFine(x2+i,y2+j)=necroticDensityFine(x2+i,y2+j)-1;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Die:
    if numel(idNec) > 0
        r2 = rand(size(idNec));
        z2 = tumorDataArray(idNec,6); % Find M nutrient level
        mMinVec = Mmin.*ones(size(z2));
        P_nec = exp(-((z2-mMinVec)/theta_nec).^2);
        go2 = (r2 < P_nec) | (z2 < mMinVec) == -1;
        idTemp2 = find(go2);
        idNecFinal = idNec(idTemp2);
        for k=1:numel(idNecFinal)
            index=idNecFinal(k);
            x=tumorDataArray(index,2);
            y=tumorDataArray(index,3);
            x1=tumorDataArray(index,4);
            y1=tumorDataArray(index,5);
            if x~=0 && y~=0
                if tdf(x,y).value > 0
                    c = tdf(x,y).value;
                    d = tdf(x,y).id(c);
                    tdf(x,y).value = c-1;
                    tdf(x,y).id(c) = 0;
                    x2 = tumorDataArray(numberOfTumorCells,2);
                    y2 = tumorDataArray(numberOfTumorCells,3);
                    f = find(tdf(x2,y2).id==numberOfTumorCells);
                    tdf(x2,y2).id(f) = d;
                    necroticDensityFine(x,y)=necroticDensityFine(x,y)+1;
                    tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)-(tumorDensityCoarse(x1,y1)>0);
                    tumorDataArray(index,1) = 0; % Mark cell as necrotic
                    temp=tumorDataArray(index,:);
                    tumorDataArray(index,:)=tumorDataArray(numberOfTumorCells,:);
                    tumorDataArray(numberOfTumorCells,:)=temp;
                    numberOfTumorCells=numberOfTumorCells-1;
                end
            end
        end
    end
    
    % Divide:
    if numel(idDiv) > 0
        r3 = rand(size(idDiv));
        z3 = tumorDataArray(idDiv,7); % Find N nutrient level
        nMinVec = Nmin.*ones(size(z3));
        P_div = 1-exp(-((z3-nMinVec)/theta_div).^2);
        go3 = r3 < P_div & z3 > nMinVec & tumorDataArray(idDiv,1) == -1;
        idTemp3 = find(go3);
        idDivFinal = idDiv(idTemp3);
        for k=1:numel(idDivFinal)
            index=idDivFinal(k);
            x=tumorDataArray(index,2);
            y=tumorDataArray(index,3);
            Tmin = tmax;
            i_min = 0; % Location with minimum tumor count
            j_min = 0;
            count = 1;
            space = zeros(8,2);
            if x~=0 && y~=0
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
                        if ~(i==0 && j==0) && tdf(x2+i,y2+j).value==0 && immuneDensityFine(x2+i,y2+j)==0
                            space(count,:) = [x2+i y2+j];
                            count = count + 1;
                        end
                    end
                end
                if numel(find(space)) > 0 % If there is a free grid space, choose a random one
                    spaceFinal = space(find(space(:,1)),:);
                    s = size(spaceFinal);
                    ri = randi([1 s(1)]);
                    chosenX = spaceFinal(ri,1);
                    chosenY = spaceFinal(ri,2);
                    tdf(chosenX,chosenY).value = 1;
                    tdf(chosenX,chosenY).id(1) = numberOfTumorCells+1;
                    x1_new=round(0.5 + h*(chosenX-0.5)); 
                    y1_new=round(0.5 + h*(chosenY-0.5));
                    tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
                    tumorDataArray(numberOfTumorCells+1,1:5)=[-1 chosenX chosenY x1_new y1_new];
                    numberOfTumorCells = numberOfTumorCells + 1;
                    if hostDensityFine(chosenX,chosenY) > 0
                        hostDensityFine(chosenX,chosenY)=0;
                        hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
                    elseif necroticDensityFine(chosenX,chosenY) > 0 % necrotic cell
                        necroticDensityFine(chosenX,chosenY)=necroticDensityFine(chosenX,chosenY)-1;
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
                            Tcount=tdf(x2+i,y2+j).value;
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
                        c = tdf(x2+i_min,y2+j_min).value;
                        tdf(x2+i_min,y2+j_min).value = c+1;
                        tdf(x2+i_min,y2+j_min).id(c+1) = numberOfTumorCells+1;
                        x1_new=round(0.5 + h*(x2+i_min-0.5));
                        y1_new=round(0.5 + h*(y2+j_min-0.5));
                        tumorDensityCoarse(x1_new,y1_new)=tumorDensityCoarse(x1_new,y1_new)+1;
                        tumorDataArray(numberOfTumorCells+1,1:5)=[-1 x2+i_min y2+j_min x1_new y1_new];
                        numberOfTumorCells = numberOfTumorCells + 1;
                        if hostDensityFine(x2+i_min,y2+j_min) > 0
                            hostDensityFine(x2+i_min,y2+j_min)=0;
                            hostDensityCoarse(x1_new,y1_new)=hostDensityCoarse(x1_new,y1_new)-1;
                        elseif necroticDensityFine(x2+i_min,y2+j_min) > 0
                            necroticDensityFine(x2+i_min,y2+j_min)=necroticDensityFine(x2+i_min,y2+j_min)-1;
                        end
                    end
                end
            end
        end
    end
    
else % Immune cells   
    % Determine which cells have tumor cells nearby:
    tCount = zeros([numberOfImmuneCells 1]);
    k = 1;
    while k <= numberOfImmuneCells
        if immuneDataArray(k,1) < 0
            x=immuneDataArray(k,2);
            y=immuneDataArray(k,3);
            T_near = 0;
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
                    T_near=T_near+tdf(x2+i,y2+j).value;
                    tCount(k) = T_near;
                end
            end
            k = k+1;
        end
    end
    idKill = find(tCount>0);
    others = find(tCount==0);
    action = randi([1 2],[numel(others) 1]);
    idMig = find(action==1);
    idDie = find(action==2);
    
    % Migrate:
    for k=1:numel(idMig)
        index = idMig(k);
        iType = immuneDataArray(index,1);
        x = immuneDataArray(index,2);
        y = immuneDataArray(index,3);
        x1 = immuneDataArray(index,4);
        y1 = immuneDataArray(index,5);
        if ~(x==0 && y==0)
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
                        immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-(immuneDensityCoarse(x1,y1)>0);
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
                        if tdf(x2+i,y2+j).value > tumorMax
                            tumorMax = tdf(x2+i,y2+j).value;
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
                immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-(immuneDensityCoarse(x1,y1)>0);
                immuneDataArray(index,2:5) = [x2+i_max y2+j_max x1_new y1_new];
            end
        end
    end
    
    % Die:
    r2 = rand(size(idDie));
    z2 = tCount(idDie);
    P_LD = 1-exp(-(thetaLD./z2).^2);
    go2 = r2 < P_LD;
    idTemp2 = find(go2);
    idDieFinal = idDie(idTemp2);
     for k=1:numel(idDieFinal)
        index = idDieFinal(k);
        iType = immuneDataArray(index,1);
        x = immuneDataArray(index,2);
        y = immuneDataArray(index,3);
        x1 = immuneDataArray(index,4);
        y1 = immuneDataArray(index,5);
        if iType == -1 % Only CTL cells die by this means.
            immuneDensityFine(x,y)=0;
            immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-(immuneDensityCoarse(x1,y1)>0);
            immuneDataArray(index,:)=immuneDataArray(numberOfImmuneCells,:);
            immuneDataArray(numberOfImmuneCells,:)=[0 0 0 0 0 0 0 0];
            numberOfImmuneCells=numberOfImmuneCells-1;
        end
     end
    
     % Kill:
     for k=1:numel(idKill)
        index = idKill(k);
        iType = immuneDataArray(index,1);
        x = immuneDataArray(index,2);
        y = immuneDataArray(index,3);
        x1 = immuneDataArray(index,4);
        y1 = immuneDataArray(index,5);
        if ~(x==0 && y==0)
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
                    if tdf(x2+i,y2+j).value > 0 && numberOfTumorCells > 0
                        chosenX = x2+i;
                        chosenY = y2+j;
                        % Count number of immune cells in neighborhood of tumor cell:
                        I_near = 0;
                        for q=-1:1;
                            if chosenX+q<1
                                chosenX2=n*m+1;
                            elseif x+q>n*m
                                chosenX2=2;
                            else
                                chosenX2=chosenX;
                            end
                            for u=-1:1
                                if chosenY+u<1
                                    chosenY2=chosenY+1;
                                elseif chosenY+u>n*m
                                    chosenY2=chosenY-1;
                                else
                                    chosenY2=chosenY;
                                end
                                I_near=I_near+immuneDensityFine(chosenX2+q,chosenY2+u);
                            end
                        end
                        % Calculate probability:
                        P_imdth = 1-exp(-(I_near)^2);
                        r = rand();
                        if r < P_imdth % Tumor cell dies                           
                            c = tdf(chosenX,chosenY).value;
                            d = tdf(x,y).id(c);
                            tdf(x,y).value = c-1;
                            tdf(x,y).id(c) = 0;
                            x2 = tumorDataArray(numberOfTumorCells,2);
                            y2 = tumorDataArray(numberOfTumorCells,3);
                            f = find(tdf(x2,y2).id==numberOfTumorCells);
                            tdf(x2,y2).id(f) = d;
                            necroticDensityFine(chosenX,chosenY)=necroticDensityFine(chosenX,chosenY)+1;
                            tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)-(tumorDensityCoarse(x1,y1)>0);
                            tumorDataArray(index,1) = 0; % Mark cell as necrotic
                            temp=tumorDataArray(index,:);
                            tumorDataArray(index,:)=tumorDataArray(numberOfTumorCells,:);
                            tumorDataArray(numberOfTumorCells,:)=temp;
                            numberOfTumorCells=numberOfTumorCells-1;
                            if iType==-2 % If NK cell
                                immuneDataArray(index,:)=[-1 x y x1 y1 0 0 killLimit]; % Replace with CTL cell
                                numberOfNKCells = numberOfNKCells-1;
                            else
                                immuneDataArray(index,8) = immuneDataArray(index,8)-1;
                                if immuneDataArray(index,8)==0 % If CTL has no kills left
                                    immuneDensityFine(x,y)=0;
                                    immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-(immuneDensityCoarse(x1,y1)>0);
                                    immuneDataArray(index,:)=immuneDataArray(numberOfImmuneCells,:);
                                    immuneDataArray(numberOfImmuneCells,:)=[0 0 0 0 0 0 0 0];
                                    numberOfImmuneCells=numberOfImmuneCells-1;
                                end
                                for o=-1:1
                                    for p=-1:1
                                        if tdf(x+o,y+p).value==0 && immuneDensityFine(x+o,y+p)==0 && ~(o==0 && p==0)
                                            r = rand();
                                            T_near = tCount(index);
                                            P_L=exp(-(thetaL/T_near)^2);
                                            if r < P_L % Recruit new CTL cell to that location
                                                x1=round(0.5 + h*(x+o-0.5)); y1=round(0.5 + h*(y+p-0.5));
                                                immuneDensityFine(x+o,y+p)=1;
                                                immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)+1;
                                                immuneDataArray(numberOfImmuneCells+1,:)=[-1 x+o y+p x1 y1 0 0 killLimit];
                                                numberOfImmuneCells = numberOfImmuneCells+1;
                                            end
                                            break
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
     end
end



end

