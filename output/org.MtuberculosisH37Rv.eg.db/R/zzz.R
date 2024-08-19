datacache <- new.env(hash=TRUE, parent=emptyenv())

org.MtuberculosisH37Rv.eg <- function() showQCData("org.MtuberculosisH37Rv.eg", datacache)
org.MtuberculosisH37Rv.eg_dbconn <- function() dbconn(datacache)
org.MtuberculosisH37Rv.eg_dbfile <- function() dbfile(datacache)
org.MtuberculosisH37Rv.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.MtuberculosisH37Rv.eg_dbInfo <- function() dbInfo(datacache)

org.MtuberculosisH37Rv.egORGANISM <- "Mycobacterium tuberculosisH37Rv"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.MtuberculosisH37Rv.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.MtuberculosisH37Rv.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.MtuberculosisH37Rv.eg_dbconn())
}

