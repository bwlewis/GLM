Generalized linear models, abridged.
===============

See our old notes on this project at: http://bwlewis.github.io/GLM.

Very revised and updated notes are in process as of October, 2018.


Here are some slides from a talk at the Cleveland R User Group: https://bwlewis.github.io/GLM/October_2018_CLERUG.html

Our current experimental reference SVD-based GLM implementation, based in turn
on Dianne O'Leary's QR implementation from 1990, can be found here:
https://github.com/bwlewis/GLM/blob/master/glm.svd.r That code is robust and
fast and replicates R's column subset selection routine by falling back to R's
default rank-revealing QR factorization in edge cases. But it's still a work in
progress.
