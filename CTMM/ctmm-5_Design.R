
###########################################################################
# Study design of ctmm methods ############################################
###########################################################################

# Install package and dependencies: ---------------------------------------

# install.packages("remotes")
remotes::install_github("ecoisilva/movedesign", force = TRUE)

# or:
install.packages(c("shiny", "fontawesome"))
shiny::runGitHub("movedesign", "ecoisilva", ref = "main")

# Run Shiny app: ----------------------------------------------------------

library(movedesign)
movedesign::run_app()
