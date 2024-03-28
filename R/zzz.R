.onLoad <- function(libname, pkgname) {
  library.dynam("d_eig", pkgname, libname, now=TRUE)
  library.dynam("locpoly", pkgname, libname, now=TRUE)
  library.dynam("testing-cfuns", pkgname, libname, now=TRUE)
  # .C("HsStart")
  invisible()
}

.onUnLoad <- function(libpath) {
  library.dynam.unload("d_eig", pkgname, libname)
  library.dynam.unload("locpoly", pkgname, libname)
  library.dynam.unload("testing-cfuns", pkgname, libname)
  invisible()
}