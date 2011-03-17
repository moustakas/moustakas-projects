function read_sc1120, linefitnodust=linefitnodust, speclinefile=speclinefile, $
  snrcuts=snrcuts, hiionly=hiionly, hiiloose=hiiloose, verbose=verbose, $
  sigcut=sigcut, _extra=extra
; jm04jun05uofa
; jm04dec02uofa - updated with S/N cuts    

; HIILOOSE - return both HII and "Unknown" galaxies
    
    if (n_elements(sigcut) eq 0L) then sigcut = 5.0
       
    path = sc1120_path(/analysis)
    speclinefile = 'sc1120_speclinefit.fits.gz'
    speclinefilenodust = 'sc1120_speclinefit_nodust.fits.gz'

    linefit = mrdfits(path+speclinefile,1,/silent)
    if arg_present(linefitnodust) then linefitnodust = mrdfits(path+speclinefilenodust,1,/silent)

    ngalaxy = n_elements(linefit)
    
    if keyword_set(snrcuts) then begin

       keep = where((linefit.oii_3727[0]/linefit.oii_3727[1] gt sigcut) and $
         (linefit.oiii_5007[0]/linefit.oiii_5007[1] gt sigcut) and $
         (linefit.h_beta[0]/linefit.h_beta[1] gt sigcut))

;      oiikeep = where(linefit.oii_3727[0]/linefit.oii_3727[1] gt sigcut,noii)
;      oiiikeep = where(linefit.oiii_5007[0]/linefit.oiii_5007[1] gt sigcut,noiii)
;      hbkeep = where(linefit.h_beta[0]/linefit.h_beta[1] gt sigcut,nhb)
;      keep = cmset_op(cmset_op(oiikeep,'AND',oiiikeep),'AND',hbkeep)

;      niikeep = where(linefit.nii_6584[0]/linefit.nii_6584[1] gt sigcut,nnii)
;      hakeep = where(linefit.h_alpha[0]/linefit.h_alpha[1] gt sigcut,nha)
;      hbkeep = where(linefit.h_beta[0]/linefit.h_beta[1] gt sigcut,nhb)
;      keep = cmset_op(cmset_op(niikeep,'AND',hakeep),'AND',hbkeep)

       nkeep = n_elements(keep)

       linefit = linefit[keep]
       if arg_present(linefitnodust) then linefitnodust = linefitnodust[keep]
       
       if keyword_set(verbose) then begin

          splog, string(nkeep,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
            ' galaxies satisfy S/N > '+string(sigcut,format='(I0)')+' criterion.'

       endif
       
       ngalaxy = nkeep
       
    endif
    
    if keyword_set(hiionly) then begin
    
       if tag_exist(linefit[0],'BPT_PURE_NII_CLASS') then begin

          hii = where(strtrim(linefit.bpt_pure_nii_class,2) eq 'HII',nhii) ; NOTE!
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

    endif
    
    if keyword_set(hiiloose) then begin
    
       if tag_exist(linefit[0],'BPT_PURE_NII_CLASS') then begin

          hiiloose = where((strtrim(linefit.bpt_pure_nii_class,2) eq 'HII') or $
            ((strtrim(linefit.bpt_class,2) eq 'Unknown') and (strtrim(linefit.bpt_pure_nii_class,2) ne 'AGN')),nhiiloose)
          if (nhiiloose ne 0L) then begin
             if keyword_set(verbose) then begin
                splog, 'Retaining '+string(nhiiloose,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
                  ' HII/UNKNOWN galaxies.'
             endif
             linefit = linefit[hiiloose]
             if arg_present(linefitnodust) then linefitnodust = linefitnodust[hiiloose]
          endif else begin
             splog, 'No HII/UNKNOWN galaxies found.'
             return, -1L
          endelse
          
       endif

    endif
    
return, linefit
end
