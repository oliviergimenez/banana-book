Note for readers

When running some of the NIMBLE code in this book, you may see the warning:

[Warning] Detected backwards indexing in ...
This is likely unintended and will likely not produce valid model code.

This happens when an individual is caught on the last capture occasion. For example, if there are T = 7 occasions and first[i] = 7, then the loop used to compute the likelihood

for (j in (first[i] + 1):T) 

becomes 8:7, which triggers the backwards indexing warning.

This is an R issue that carries over when the code is compiled in C++. In principle, one would wrap the loop in a conditional such as if (first[i] < T) {...}, but such statements are not (yet) supported in NIMBLE. Another option would be to remove individuals first caught on the last occasion, though I prefer not to alter the data.

The important point is that this warning has no impact on parameter estimates. I checked this with simulations and by comparison with maximum likelihood results from standard software (e.g. MARK, E-SURGE). Individuals first caught on the last occasion do not contribute to detection or transition parameters, since their capture history ends immediately after capture.

In short: you can safely ignore this warning. If it bothers you, simply add

nimbleOptions(verbose = FALSE)

at the start of your code.

Olivier, August 22, 2025

