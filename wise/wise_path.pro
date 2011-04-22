function wise_path, ppxf=ppxf, wisesfrs=wisesfrs, catalogs=catalogs, $
  dimage=dimage
; jm11apr13ucsd
    
    wisepath = getenv('IM_RESEARCH_DIR')+'/projects/wise/'
    if keyword_set(wisesfrs) then wisepath = wisepath+'wisesfrs/'
    if keyword_set(catalogs) then wisepath = wisepath+'catalogs/'
    if keyword_set(dimage) then wisepath = getenv('DIMAGE_DIR')+'/data/atlas/'
    if keyword_set(ppxf) then wisepath = wisepath+'ppxf/'

return, wisepath
end
