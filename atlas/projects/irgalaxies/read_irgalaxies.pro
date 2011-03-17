function read_irgalaxies, linefitnodust=linefitnodust, speclinefile=speclinefile, $
  snrcuts=snrcuts, hiionly=hiionly, verbose=verbose, niikeep=niikeep, hakeep=hakeep, $
  hbkeep=hbkeep, oiikeep=oiikeep, oiiikeep=oiiikeep, rasort=rasort, _extra=extra, $
  templates_04=templates_04
; jm05aug12uofa
; TEMPLATES_04 - read the fitting results from just 4 templates    
    
    path = atlas_path(/projects)+'irgalaxies/specfit/'
    if keyword_set(templates_04) then begin
       speclinefile = 'ir_04_speclinefit.fits.gz'
       speclinefilenodust = 'ir_04_speclinefit_nodust.fits.gz'
    endif else begin
       speclinefile = 'ir_integrated_speclinefit.fits.gz'
       speclinefilenodust = 'ir_integrated_speclinefit_nodust.fits.gz'
    endelse

    linefit = mrdfits(path+speclinefile,1,/silent)
    if arg_present(linefitnodust) then linefitnodust = mrdfits(path+speclinefilenodust,1,/silent)

    ngalaxy = n_elements(linefit)

    cut = 3.0
    
    niikeep = where(linefit.nii_6584[0]/linefit.nii_6584[1] ge cut,nnii,comp=nonii)
    hakeep = where(linefit.h_alpha[0]/linefit.h_alpha[1] ge cut,nha,comp=noha)
    hbkeep = where(linefit.h_beta[0]/linefit.h_beta[1] ge cut,nhb,comp=nohb)
    oiikeep = where(linefit.oii_3727[0]/linefit.oii_3727[1] ge cut,noii,comp=nooii)
    oiiikeep = where(linefit.oiii_5007[0]/linefit.oiii_5007[1] ge cut,noiii,comp=nooiii)

    if keyword_set(snrcuts) then begin

       keep = cmset_op(hakeep,'AND',hbkeep)
;      keep = cmset_op(cmset_op(cmset_op(niikeep,'AND',hakeep),'AND',hbkeep),'AND',oiikeep)
       nkeep = n_elements(keep)

       linefit = linefit[keep]
       if arg_present(linefitnodust) then linefitnodust = linefitnodust[keep]
       
       if keyword_set(verbose) then begin
          splog, string(nkeep,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
            ' total galaxies satisfy S/N > 3 criterion.'
       endif
       
       ngalaxy = nkeep
       
    endif
    
    if keyword_set(hiionly) then begin

       hii = where((strtrim(linefit.bpt_class,2) eq 'HII'),nhii) ; NOTE!
       unknown = where((strtrim(linefit.bpt_class,2) eq 'Unknown'),nunknown)
       
       if (nhii ne 0L) then begin
          if keyword_set(verbose) then begin
             splog, 'Retaining '+string(nhii,format='(I0)')+'/'+$
               string(ngalaxy,format='(I0)')+' HII galaxies.'
          endif
          linefit = linefit[hii]
          if arg_present(linefitnodust) then linefitnodust = linefitnodust[hii]
       endif else begin
          splog, 'No HII galaxies found.'
          return, -1L
       endelse
       
    endif

    if keyword_set(rasort) then begin
       srt = sort(linefit.ra)
       linefit = linefit[srt]
       if arg_present(linefitnodust) then linefitnodust = linefitnodust[srt]
    endif
    
return, linefit
end
