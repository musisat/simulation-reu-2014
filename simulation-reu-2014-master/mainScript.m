% main script 
global tdf necroticDensityFine maxNumberOfTumorCells
global numberOfImmuneCells immuneDataArray nutrientM n m
global nutrientN numberOfTumorCells tumorDataArray

initializeVariables;
I0 = initializeCells(10);
iters = 3;
tumorHistory = zeros(1,iters);
for generation = 1:iters
    if numberOfTumorCells > maxNumberOfTumorCells
        error('Maximum number of cells exceeded.')
    elseif numberOfTumorCells == 0
        error('No more tumor cells left');
    end
    nutrientM=calculateNutrientM;
    nutrientN=calculateNutrientN;
    for k=1:numberOfTumorCells
        x = tumorDataArray(k,2);
        y = tumorDataArray(k,3);
        x1 = tumorDataArray(k,4);
        y1 = tumorDataArray(k,5);
        tumorDataArray(k,6) = nutrientM(x1,y1);
        tumorDataArray(k,7) = nutrientN(x1,y1);
    end
    for k=1:numberOfImmuneCells
        x1 = immuneDataArray(k,4);
        y1 = immuneDataArray(k,5);
        immuneDataArray(k,6)=nutrientM(x1,y1);
        immuneDataArray(k,7)=nutrientN(x1,y1);
    end
    cellActionVec(1);
    cellActionVec(2);
    produceNK(I0);
    tumorHistory(generation) = numberOfTumorCells;
end

figure(1);
subplot(2,2,1);
tumorDensityFine = zeros(n*m);
for x=1:n*m
    for y=1:n*m
        tumorDensityFine(x,y)=tdf(x,y).value;
    end
end
imagesc(tumorDensityFine); hold on;
for k=1:numberOfImmuneCells
   pt = immuneDataArray(k,2:3);
   iType = immuneDataArray(k,1);
   if iType==-1
     plot(pt(1),pt(2),'m*','MarkerSize',10);
   else
     plot(pt(1),pt(2),'r*','MarkerSize',10);
   end
end
title('Tumor Density with Immune Cells (stars)'); 
hold off;
subplot(2,2,2);
imagesc(nutrientN);
subplot(2,2,3);
imagesc(necroticDensityFine);
subplot(2,2,4);
plot(tumorHistory);


