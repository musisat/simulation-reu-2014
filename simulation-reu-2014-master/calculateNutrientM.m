function M = calculateNutrientM() 
   global lambdaTM L source n tumorDensityCoarse immuneDensityCoarse
   global hostDensityCoarse lambdaH

   %beta=lambdaTM*alpha^2*tumorDensityCoarse + alpha^2*immuneDensityCoarse + alpha^2*hostDensityCoarse; 
   beta=lambdaTM*tumorDensityCoarse + lambdaH*immuneDensityCoarse + lambdaH*hostDensityCoarse; 
   beta=beta(:); 
   A = L + spdiags(beta,0,n^2,n^2); 
   M = A\source;
   M = reshape(M,n,n); 
   
end

