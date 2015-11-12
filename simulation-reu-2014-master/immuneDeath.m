function immuneDeath(x,y,x1,y1,index)

global immuneDataArray immuneDensityFine immuneDensityCoarse 
global numberOfImmuneCells

immuneDensityFine(x,y)=0;
immuneDensityCoarse(x1,y1)=immuneDensityCoarse(x1,y1)-1;
immuneDataArray(index,:)=immuneDataArray(numberOfImmuneCells,:);
immuneDataArray(numberOfImmuneCells,:)=[0 0 0 0 0 0 0 0];
numberOfImmuneCells=numberOfImmuneCells-1;

end

