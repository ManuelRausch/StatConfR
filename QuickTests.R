################################################################################
######### ANALYZE R-PACKAGE STATCONFR ##########################################
################################################################################

# Manuel Rausch, 15.04.2023
# We tested the full functionality of the package with the script "TestScript.R".
# However, running that script takes a considerable amount of time.
# QuickTests.R" is a script for those who need to test whether the package is working in principle, but
# do not have the time to run "TestScript.R".

# Install the development version of the package

devtools::install_github("ManuelRausch/StatConfR")

devtools::check(StatConfR)

cran_checks <- rhub::check_for_cran()
cran_checks$cran_summary()


