# Create a compressed version for the dataset 'DATA'
tools::checkRdaFiles('~/Dropbox/Repos/fuse/data/DATA.rda')
tools::resaveRdaFiles(paths = '~/Dropbox/Repos/fuse/data/DATA.rda',
                      compress = 'bzip2')

# Create a compressed version for the dataset 'modlist'
tools::checkRdaFiles('~/Dropbox/Repos/fuse/data/modlist.rda')
tools::resaveRdaFiles(paths = '~/Dropbox/Repos/fuse/data/modlist.rda',
                      compress = 'bzip2')

# Create a compressed version for the dataset 'modliststring'
tools::checkRdaFiles('~/Dropbox/Repos/fuse/data/modliststring.rda')
tools::resaveRdaFiles(paths = '~/Dropbox/Repos/fuse/data/modliststring.rda',
                      compress = 'gzip')

# Create a compressed version for the dataset 'parameters'
tools::checkRdaFiles('~/Dropbox/Repos/fuse/data/parameters.rda')
tools::resaveRdaFiles(paths = '~/Dropbox/Repos/fuse/data/parameters.rda',
                      compress = 'gzip')

# Create the Appveyor config file for continuous integration on Windows
devtools::use_appveyor()
# move the newly created appveyor.yml to the root directory and modify it

# Create the travis config file for continuous integration on Linux-OSX
devtools::use_travis()

# Generate a template for a README.Rmd
devtools::use_readme_rmd()

# Generate a template for a Code of Conduct
devtools::use_code_of_conduct()

# Check spelling mistakes
devtools::spell_check('fuse',
                      ignore = c('metres', 'catalogue', 'DEFRA', 'EMEP',
                                 'EPSG', 'WGS'))

# Run R CMD check
devtools::check('fuse')
# The above will also run the unit tests using testthat
# devtools::test('fuse')
