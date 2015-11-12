function cellAction(x,y,x1,y1,cellType,index)
%CELLACTION will decide which action each cell will take and update the
%   grid each time an action is taken. Inputs and outputs are mostly global
%   variables (for now). cellType = 1 for immune cells and 2 for tumor
%   cells. For now, assume that a singular action can occur at each grid
%   location at each time step, and the action depends on the most
%   prominent cell type that occupies that location.

global immuneDataArray tumorDensityFine thetaLD killLimit tumorDataArray
global immuneDensityFine numberOfNKCells theta_nec nutrientN nutrientM
global numberOfTumorCells theta_div theta_mig localM localN n m Nmin Mmin

localM = nutrientM(x1,y1);
localN = nutrientN(x1,y1);
r = rand(); % Random number for coinflips

if cellType==1 % For immune cell
% Determine immune cell type:
iType = immuneDataArray(index,1);
T_near = 0; % Tumor cells in neighborhood
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
        T_near=T_near+tumorDensityFine(x2+i,y2+j);
    end
end
if T_near==0 % If no tumor cells around, either die or migrate
    if iType==-2 % If NK cell
        immuneAction = 2; % Migrate is only option
    else
        immuneAction = randi([1 2]); % CTL cell can either die or migrate
    end
    if immuneAction==1 % Cell marked for death
        P_LD = 1-exp(-(thetaLD/T_near)^2);       
        if r < P_LD
            immuneDeath(x,y,x1,y1,index);
        end
    else % Cell marked for migration
        immuneMigrate(x,y,x1,y1,index,iType)
    end   
else % Immune cell attempts kill
    % Choose nearby tumor cell:
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
        if tumorDensityFine(x2+i,y2+j) > 0
            chosenX = x2+i;
            chosenY = y2+j;
            done = 1;
        end
    end
    % Count number of immune cells in neighborhood of tumor cell:
    I_near = 0;
    for i=-1:1
        if chosenX+i<1
            chosenX2=n*m+1;
        elseif x+i>n*m
            chosenX2=2;
        else
            chosenX2=chosenX;
        end
        for j=-1:1
            if chosenY+j<1
                chosenY2=chosenY+1;
            elseif chosenY+j>n*m
                chosenY2=chosenY-1;
            else
                chosenY2=chosenY;
            end
            I_near=I_near+immuneDensityFine(chosenX2+i,chosenY2+j);
        end
    end
    % Calculate probability:
    P_imdth = 1-exp(-(I_near)^2);
    if r < P_imdth % Tumor cell dies
        %immdth = 1
        indices = find(tumorDataArray(:,1)==-1 & tumorDataArray(:,2)==chosenX & tumorDataArray(:,3)==chosenY);
        if numel(indices) == tumorDensityFine(chosenX,chosenY)
            g = randi([1 tumorDensityFine(chosenX,chosenY)]);
            tIndex = indices(g);
            tumorNecrosis(x,y,x1,y1,tIndex);
            if iType==-2 % If NK cell
                immuneDataArray(index,:)=[-1 x y x1 y1 0 0 killLimit]; % Replace with CTL cell
                numberOfNKCells = numberOfNKCells-1;
            else
                immuneDataArray(index,8) = immuneDataArray(index,8)-1;
                if immuneDataArray(index,8)==0 % If CTL has no kills left
                    immuneDeath(x,y,x1,y1,index);
                end
                induceCTL(x,y); % When CTL makes kill, induce more CTLs
            end
        end
    end
end

else % For tumor cells:
tAction = randi([1 3]); % 1=die, 2=divide, 3=move
if tAction==1
    %P_nec = exp(-(localM/(numberOfTumorCells*theta_nec))^2);
    P_nec = exp(-((localM-Mmin)/theta_nec)^2);
    if r < P_nec
        tumorNecrosis(x,y,x1,y1,index);
    end
elseif tAction==2
    %P_div = 1-exp(-(localN/(numberOfTumorCells*theta_div))^2);
    P_div = 1-exp(-((localN-Nmin)/theta_div)^2);
    if r < P_div
        tumorDivide(x,y);
    end
else
    P_mig = 1-exp(-numberOfTumorCells*(localM/theta_mig)^2);
    if r < P_mig
        tumorMigrate(x,y,x1,y1,index);
    end
end
end

% z = tumorDataArray(1:numberOfTumorCells,8)
% probDiv = ...... (vectorized)i_min = 0; % Location with minimum tumor count
% r = rand(1:numberOfTumorCells,1)
% divide = r < probDiv
% id = find(divide)
% for k = id
   

end

