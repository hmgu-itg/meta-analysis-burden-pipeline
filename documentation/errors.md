# Mismatching .score and .var file used
This error is raised when the `.var.{i}` file is larger than expected during the `split_scorevar` rule. 

`.score.{i}` file and the corresponding `.var.{i}` file is read to recreate all the covariance matrix. The covariance matrix for the groups in the `.score.{i}` is recreated sequentially based on the values read from the `.var.{i}` file. The number of values fetched from the `.var.{i}` is determined by the number of variants for the group (in the `.score.{i}` table). Even after recreating all covariance matrix and there's still values left in the `.var.{i}` file, it indicates that the input score and var files are not the matching pair of each other, and so the error is raised.  

If you encounter this issue, you should either rerun the single-cohort burden testing pipeline to regenerate those files, or exclude them from the analysis.

# Error in group.V[lower.tri(group.V, diag = TRUE)] <- readBin...

`Error in group.V[lower.tri(group.V, diag = TRUE)] <- readBin(con, what = "numeric", : (converted from warning) number of items to replace is not a multiple of replacement length`

This error is raised when the `.var.{i}` file is smaller than expected during the `split_scorevar` rule.

As explained in the `Mismatching .score and .var file used` error, we have an expectation of what the size of `.var.{i}` file would be based on the contents of the corresponding `.score.{i}` file. Encountering this error during recreating the covariance matrix from the values in the `.var.{i}` file indicates that there's no more values to read from the `.var.{i}` and so the covariance matrix cannot be recreated. 

Again, if you encounter this issue, you should either rerun the single-cohort burden testing pipeline to regenerate those files, or exclude them from the analysis.