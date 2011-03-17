function sfh_get_habs, ppxf
; jm10may06ucsd - compute (Hg+Hd)/2 under a variety of conditions

    common get_habs, lick

    ngal = n_elements(ppxf)
    
; read the Lick absorption-line table    
    if (n_elements(lick) eq 0) then begin
       indexfile = getenv('IMPRO_DIR')+'/etc/indexlist.dat'
       readcol, indexfile, licknames, w1, w2, wb1, wb2, $
         wr1, wr2, units, format='A,D,D,D,D,D,D,A', /silent, $
         comment='#'
       lick = replicate({index: '', w1: 0.0, w2: 0.0, wb1: 0.0, $
         wb2: 0.0, wr1: 0.0, wr2: 0.0, units: ''},n_elements(w1))
       lick.index = licknames & lick.w1 = w1 & lick.w2 = w2
       lick.wb1 = wb1 & lick.wb2 = wb2 & lick.wr1 = wr1
       lick.wr2 = wr2 & lick.units = units
    endif

; read the list of skylines
    readcol, getenv('IMPRO_DIR')+'/etc/skylines.dat', $
      skywave, format='A', comment='#'
    nsky = n_elements(skywave)
    
; check for any skyline falling in any of the Hd,Hg index windows  
    z1 = 1.0+ppxf.z
    for ii = 0, nsky-1 do begin
       hg = where(strmatch(lick.index,'*lick_hg_a*',/fold))
       hgcrap = where($
         ((skywave[ii] ge lick[hg].w1*z1)  and (skywave[ii] le lick[hg].w2*z1)) or $
         ((skywave[ii] ge lick[hg].wb1*z1) and (skywave[ii] le lick[hg].wb2*z1)) or $
         ((skywave[ii] ge lick[hg].wr1*z1) and (skywave[ii] le lick[hg].wr2*z1)),nhgcrap)
       hd = where(strmatch(lick.index,'*lick_hd_a*',/fold))
       hdcrap = where($
         ((skywave[ii] ge lick[hd].w1*z1)  and (skywave[ii] le lick[hd].w2*z1)) or $
         ((skywave[ii] ge lick[hd].wb1*z1) and (skywave[ii] le lick[hd].wb2*z1)) or $
         ((skywave[ii] ge lick[hd].wr1*z1) and (skywave[ii] le lick[hd].wr2*z1)),nhdcrap)
       if (nhgcrap gt 0) or (nhdcrap gt 0) then begin
          splog, 'Skyline '+string(skywave[ii],format='(I0)')
          splog, '  Hg: '+strtrim(nhgcrap,2)+'/'+strtrim(ngal,2)
          splog, '  Hd: '+strtrim(nhdcrap,2)+'/'+strtrim(ngal,2)
       endif
    endfor
    
stop    
      
; default: average of the two indices    
    avg = (ppxf.lick_hg_a_cor[0]+ppxf.lick_hd_a_cor[0])/2.0


stop    
return, habs
end
