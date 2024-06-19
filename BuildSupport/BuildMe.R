
# Pulls from repo to get latest changes
system("git pull --ff-only")


if(FALSE){
  
  # reverses local commit (but changes will still be in the files)
  system("git reset HEAD^") # 
  
  # stashes the local file changes
  system("git stash")
  
  # pulls changes from remote git and merges with local
  system("git pull --ff-only")
  
  # returns your stashed local changes
  system("git stash pop")
  
}


# project-specific tasks --------------------------------------------------
# clear environment
rm(list = ls()); gc()
# .rs.restartR()
# manually update files ---------------------------------------------------

file.edit('DESCRIPTION')
file.edit('NEWS')

file.edit(file.path(getwd(), "BuildSupport", "source_search.R"))

# update and develop ------------------------------------------------------
devtools::document(); devtools::load_all()

# publication tools
rmarkdown::render("README.Rmd", 
                  output_format = "md_document",
                  output_file="README.md")


devtools::test(stop_on_failure = TRUE)

# lifecycle::last_lifecycle_warnings()

devtools::run_examples()

devtools::check_man()

# check -------------------------------------------------------------------

# check complete package
devtools::check(vignettes = TRUE )
# devtools::check(vignettes = FALSE )

# install
devtools::install(upgrade = "never")


file.edit("README.md")

pkgdown::build_site()
# usethis::use_apache_license()
# buildignore -------------------------------------------------------------

usethis::use_build_ignore(c(
  "BuildSupport"
))

