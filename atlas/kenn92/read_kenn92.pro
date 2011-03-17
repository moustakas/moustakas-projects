function read_kenn92, linefitnodust=linefitnodust, speclinefile=speclinefile, $
  snrcuts=snrcuts, hiionly=hiionly, verbose=verbose, niikeep=niikeep, hakeep=hakeep, $
  hbkeep=hbkeep, oiikeep=oiikeep, oiiikeep=oiiikeep, _extra=extra
; jm05aug03uofa
    
    path = kenn92_path(/specfit)
    speclinefile = 'kenn92_speclinefit.fits.gz'
    speclinefilenodust = 'kenn92_speclinefit_nodust.fits.gz'

    linefit = mrdfits(path+speclinefile,1,/silent)
    if arg_present(linefitnodust) then linefitnodust = mrdfits(path+speclinefilenodust,1,/silent)

    ngalaxy = n_elements(linefit)
    
    cut = 3.0
    
    niikeep = where(linefit.nii_6584[0]/linefit.nii_6584[1] gt cut,nnii,comp=nonii)
    hakeep = where(linefit.h_alpha[0]/linefit.h_alpha[1] gt cut,nha,comp=noha)
    hbkeep = where(linefit.h_beta[0]/linefit.h_beta[1] gt cut,nhb,comp=nohb)
    oiikeep = where(linefit.oii_3727[0]/linefit.oii_3727[1] ge cut,noii,comp=nooii)
    oiiikeep = where(linefit.oiii_5007[0]/linefit.oiii_5007[1] ge cut,noiii,comp=nooiii)

    if keyword_set(snrcuts) then begin

       keep = cmset_op(hakeep,'AND',hbkeep)
       nkeep = n_elements(keep)

       linefit = linefit[keep]
       if arg_present(linefitnodust) then linefitnodust = linefitnodust[keep]
       
       if keyword_set(verbose) then begin

;         splog, string(nnii,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' galaxies satisfy [N II] > 3.'
;         splog, string(nha,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' galaxies satisfy Ha > 3.'
;         splog, string(nhb,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' galaxies satisfy Hb > 3.'

          splog, string(nkeep,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' galaxies satisfy S/N > 3 criterion.'

       endif
       
       ngalaxy = nkeep
       
    endif
    
    if keyword_set(hiionly) then begin
    
       hii = where(strtrim(linefit.bpt_class,2) eq 'HII',nhii)

       if (nhii ne 0L) then begin
          if keyword_set(verbose) then begin
             splog, 'Retaining '+string(nhii,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
               ' HII galaxies.'
          endif
          linefit = linefit[hii]
          if arg_present(linefitnodust) then linefitnodust = linefitnodust[hii]
       endif else begin
          splog, 'No HII galaxies found.'
          return, -1L
       endelse
          
    endif 
    
return, linefit
end
