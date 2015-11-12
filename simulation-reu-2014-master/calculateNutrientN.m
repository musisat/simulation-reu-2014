function N = calculateNutrientN()
   global lambdaTN L source n tumorDensityCoarse 
   global hostDensityCoarse immuneDensityCoarse lambdaH

   %beta=lambdaTN*alpha^2*tumorDensityCoarse + alpha^2*immuneDensityCoarse + alpha^2*hostDensityCoarse; 
   beta=lambdaTN*tumorDensityCoarse + lambdaH*immuneDensityCoarse + lambdaH*hostDensityCoarse; 
   beta=beta(:); 
   A = L + spdiags(beta,0,n^2,n^2); 
   N = A\source;
   N = reshape(N,n,n); 
   
end



