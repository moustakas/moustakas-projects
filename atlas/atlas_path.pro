function atlas_path, atlas1d=atlas1d, atlas2d=atlas2d, ascii=ascii, web=web, $
  dataweb=dataweb, literature=literature, papers=papers, projects=projects, $
  analysis=analysis, templates1d=templates1d, templates2d=templates2d, $
  ned=ned, dss=dss, textables=textables, plots=plots, turner=turner, kenn92=kenn92, $
  iue=iue, specfit=specfit, goldmine=goldmine, atlasweb=atlasweb, $
  spec2datlas=spec2datlas, spec2dturner=spec2dturner, spectral_atlas=spectral_atlas, $
  ppxf=ppxf
; jm02apr15uofa
; jm09dec17uofa - added PPXF support

    home = getenv('RESEARCHPATH')+'/projects/atlas'
    
    path = home+'/'
    if keyword_set(spectral_atlas) then path = home+'/spectral_atlas/'
    if keyword_set(atlas1d) then path = home+'/atlas1d/'
    if keyword_set(ascii) then path = home+'/ascii/'
    if keyword_set(web) then path = home+'/public_html/spectral_atlas/'
;   if keyword_set(web) then path = getenv('BASEPATH')+'/ay/public_html/research/spectral_atlas/'
    if keyword_set(dataweb) then path = '/d1/ioannis/spec2datlas/public_html/redux/'
;   if keyword_set(dataweb) then path = home+'/public_html/atlas/redux/'
    if keyword_set(literature) then path = home+'/literature/'
    if keyword_set(papers) then path = home+'/papers/'
    if keyword_set(projects) then path = home+'/projects/'
    if keyword_set(analysis) then path = home+'/analysis/'
    if keyword_set(templates1d) then path = home+'/projects/templates/spec1d/'
    if keyword_set(templates2d) then path = home+'/projects/templates/spec2d/'
    if keyword_set(ned) then path = home+'/analysis/ned/'
    if keyword_set(atlasweb) then path = home+'/webatlas/'
    if keyword_set(dss) then path = home+'/DSS/'
    if keyword_set(textables) then path = home+'/papers/textables/'
    if keyword_set(plots) then path = home+'/papers/plots/'
    if keyword_set(turner) then path = home+'/projects/turner/'
    if keyword_set(kenn92) then path = home+'/projects/kenn92/'
    if keyword_set(iue) then path = home+'/projects/iue/'
    if keyword_set(specfit) then path = home+'/specfit/'
    if keyword_set(goldmine) then path = home+'/projects/goldmine/'
    if keyword_set(spec2dturner) then path = home+'/projects/turner/spec2dturner/'

    if keyword_set(spec2datlas) then begin
       if file_test('/global/bias1/ioannis/spectral_atlas/raw/',/directory) then begin
          path = '/global/bias1/ioannis/spectral_atlas/raw/' 
       endif else begin
          if file_test('/Volumes/WDexternal/data/spectral_atlas/raw/',/directory) then $
            path = '/Volumes/WDexternal/data/spectral_atlas/raw/' else message, 'External drive not mounted!'
       endelse
    endif
    
    if keyword_set(atlas2d) then begin
       if file_test('/global/bias1/ioannis/spectral_atlas/spec2d/',/directory) then begin
          path = '/global/bias1/ioannis/spectral_atlas/spec2d/' 
       endif else begin
          if file_test('/Volumes/WDexternal/data/spectral_atlas/spec2d/',/directory) then $
            path = '/Volumes/WDexternal/data/spectral_atlas/spec2d/' else message, 'External drive not mounted!'
       endelse
    endif

; PPXF
    if keyword_set(ppxf) then path = getenv('RESEARCHPATH')+'/data/atlas/ppxf/'
    
return, path
end
