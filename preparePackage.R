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

# Run unit tests using testthat
devtools::test('fuse')

# Run R CMD check or devtools::check()
devtools::check('fuse')

# Create the Appveyor config file for continuous integration on Windows
devtools::use_appveyor()
# move the newly created appveyor.yml to the root directory and modify it
