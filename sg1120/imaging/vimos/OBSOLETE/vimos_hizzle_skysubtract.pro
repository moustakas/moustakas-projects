pro vimos_hizzle_skysubtract, dec03=dec03, feb06=feb06, sg1120=sg1120, standards=standards
; jm07jan16nyu - sky-subtract everything

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    sg1120_path = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'
    standards_path = vimos_path(dec03=dec03,feb06=feb06)+'standards/'

    bandpass = ['B','V','R']
    quadrant = ['Q1','Q2','Q3','Q4']

    if keyword_set(sg1120) then begin
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin
             splog, '###########################################################################'
             splog, 'SG1120 '+bandpass[ib]+'-band, '+quadrant[iq]
             splog, '###########################################################################'
             imagelist = file_search(sg1120_path+'ra.sg1120*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             hizzle_skysubtract, imagelist, outdir=sg1120_path, apradius=apradius, nsamp=nsamp, $
               maskwidth=maskwidth, sigmacut=sigmacut
          endfor
       endfor
    endif
       
    if keyword_set(standards) then begin
       for ib = 0L, n_elements(bandpass)-1L do begin
          for iq = 0L, n_elements(quadrant)-1L do begin
             splog, '###########################################################################'
             splog, 'Standards '+bandpass[ib]+'-band, '+quadrant[iq]
             splog, '###########################################################################'
             imagelist = file_search(standards_path+'ra.[PG,SA]*_'+bandpass[ib]+'_'+quadrant[iq]+'.fits')
             hizzle_skysubtract, imagelist, outdir=standards_path, apradius=apradius, nsamp=nsamp, $
               maskwidth=maskwidth, sigmacut=sigmacut
          endfor
       endfor
    endif
       
return
end
