pro get_bcg_radec

;   ff = file_search('13may20/*f160w_drz*.txt',count=nn)

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    ncc = n_elements(clash)
    for ic = 0, ncc-1 do begin
;      splog, 'Cluster '+clash[ic].shortname
       if clash[ic].shortname eq 'a2261' then begin
          mosaicpath = getenv('CLASH_ARCHIVE')+'/'+clash[ic].dirname+'/HST/images/'+$
            'mosaicdrizzle_image_pipeline/scale_30mas/'
       endif else begin
          mosaicpath = getenv('CLASH_ARCHIVE')+'/'+clash[ic].dirname+'/HST/images/'+$
            'mosaicdrizzle_image_pipeline/scale_65mas/'
       endelse          
       
       ff = file_search('13may20/'+clash[ic].shortname+'_*f160w*_drz*.txt',count=nn)
;         if clash[ic].shortname eq 'a2261' then stop
       if nn eq 1 then begin
          txt = djs_readilines(ff,indx=2)
          coord = strsplit((strsplit(txt,'(',/extract))[1],'),',/extract)
          xcen = float(coord[1]) ; pixels
          ycen = float(coord[0])

          imfile = file_search(mosaicpath+clash[ic].shortname+$
            '_mosaic_???mas_*_f160w*drz_????????.fits*',count=nn)

;         splog, mosaicpath+clash[ic].shortname+'_mosaic_065mas_*_f160w_drz_????????.fits*'
          if nn ne 1 then stop
;         splog, imfile

          hdr = headfits(imfile)
          extast, hdr, astr
          xy2ad, xcen, ycen, astr, ra, dec
;         splog, clash[ic].shortname, ra, dec
       endif else splog, 'No file for cluster '+clash[ic].shortname
    endfor
       
return
end
    
