# Compute the Fractions Skill Score for ensembles

`test_fss.R` is a template for how to compute the ensemble fractions skill score. 

You need to set a number of paths and other information for both the forecast and the observations files. You also need to explicitly set 
the verification domain. In `test_fss.R` this is done by reading a single forecast file and using `subgrid()` from the {_meteogrid_} package 
to take a cut-out of the domain. 

The script calls `ens_read_and_fbs()`, which computes the unaggregated Fractions Brier Score (fbs) and the reference Fractions Brier Score (fbs_ref). 
The function is called for each forecast date and lead time by using `expand.grid()` to get every combination and `map2_dfr()` from the {_purrr_} 
package to call the function for each combination. This can be changed, though it is best to at least treat each forecast date separately. 

From the output of `ens_read_and_fbs()` the Fractions Skill Score can be calculated by grouping the data (in the example by the forecast model, 
lead time, neighbourhood length, and threshold) and computing 1 - (sum(fbs) / sum(fbs_ref)). 

A couple of methods for plotting the output are also included. 

Note that the {_here_} package and the development versions of the {_harpIO_} and {_harpVis_} packages are needed, and can be installed with:

```{r}
install.packages("here")
remotes::install_github("andrew-met/harpIO", "develop")
remotes::install_github("andrew-met/harpVis", "gg_geofield")
```
