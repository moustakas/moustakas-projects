pro kenn92_view_webatlas, object, webatlas=webatlas
; jm04jan19uofa
    
; initialize path names
    
    dsspath = atlas_path(/kenn92)+'webkenn92/DSS/' ; image data path
    webpath = atlas_path(/web)+'kenn92/'           ; web path
    specfitpath = atlas_path(/kenn92)+'analysis/'

; restore all the fitting results

    atlas = read_kenn92(/silent)

    view_webatlas, atlas, object, specfitpath=specfitpath, dsspath=dsspath, $
      webpath=webpath, webatlas=webatlas, root='kenn92'

return
end
