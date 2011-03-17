function read_literature_halpha_catalog
; jm05aug17uofa
; read Kennicutt's average literature H-alpha catalog

    catpath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='halpha')
    catfile = 'literature_halpha_catalog_05aug17.sex'

    cat = rsex(catpath+catfile)

return, cat
end
