######################################################################
#
# zzz.R
#
# Written by Christopher DuBois
#
# Last Modified 11/1/2010
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# .First.lib is run when the package is loaded with library(mpmm)
#
######################################################################

.First.lib <- function(lib, pkg){
   library.dynam("mpmm", pkg, lib)
   cat(paste("copyright (c) 2010, Christopher DuBois,",
             "University of California-Irvine\n",sep=""))
   cat('Type help(package="mpmm") to get started.\n')
}
