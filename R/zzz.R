.onLoad <- function(libname, pkgname) {
  path <- system.file("libs", package = pkgname)
  dyn.load(file.path(path, "d_eig.dll"))
  dyn.load(file.path(path, "locpoly.dll"))
  dyn.load(file.path(path, "testing-cfuns.dll"))
}
