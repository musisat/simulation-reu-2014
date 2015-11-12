function tumorNecrosis(x,y,x1,y1,id)

global tumorDataArray tdf tumorDensityCoarse numberOfTumorCells
global necroticDensityFine

c = tdf(x,y).value;
tdf(x,y).value=c-1;
index = find(tdf(x,y).id) == id;
tdf(x,y).id(index) = 0;
necroticDensityFine(x,y)=necroticDensityFine(x,y)+1;
tumorDensityCoarse(x1,y1)=tumorDensityCoarse(x1,y1)-1;
%tumorDataArray(tIndex,1:5) = [0 0 0 0 0]; % Mark cell as necrotic
%temp=tumorDataArray(tIndex,:);
%tumorDataArray(tIndex,:)=tumorDataArray(numberOfTumorCells,:);
%tumorDataArray(numberOfTumorCells,:)=temp;
numberOfTumorCells=numberOfTumorCells-1;

end

