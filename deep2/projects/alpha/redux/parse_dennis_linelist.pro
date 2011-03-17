pro parse_dennis_linelist
; jm08jun07nyu - put Dennis' ThAr linelist into xidl format

    inpath = getenv('RESEARCHPATH')+'/data/deep2_alpha/linelists/'
    outpath = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/'

; read the parent line-list to extrapolate blueward and redward
; of Dennis' line-list    
    readcol, outpath+'mike_thar.lst', xwave, xflag, xjunk, $
      format='D,A,A', /silent
    readcol, inpath+'thar_08jun10.lst', wave, flag, format='D,A', /silent

;   pre = where(xwave lt min(wave),npre)
;   wave = [xwave[pre],wave]
;   flag = [replicate('M',npre),flag]

    pre = where(xwave lt min(wave),npre)
    post = where(xwave gt max(wave),npost)
    wave = [xwave[pre],wave,xwave[post]]
    flag = [replicate('M',npre),flag,replicate('M',npost)]
    
    nall = n_elements(wave)
    these = where(flag eq 'M',nthese)

; everything    
    openw, lun, outpath+'mike_thar_alpha_custom.lst', /get_lun
    for ii = 0L, nall-1L do printf, lun, wave[ii], '1', 'ThAr', $
      format='(F9.4,1x,A1,1x,A4)'
    free_lun, lun
    
; just the murphy line-list
    openw, lun, outpath+'mike_thar_alpha_murphy.lst', /get_lun
    for ii = 0L, nthese-1L do printf, lun, wave[these[ii]], $
      '1', 'ThAr', format='(F9.4,1x,A1,1x,A4)'
    free_lun, lun

; finally, generate a line-list that has had a global shift applied,
; with no distinction between Murphy and non-Murphy lines

    readcol, inpath+'thar_08jun09_shift.lst', wave, format='D', /silent

    pre = where(xwave lt min(wave),npre)
    post = where(xwave gt max(wave),npost)
    wave = [xwave[pre],wave,xwave[post]]
    nall = n_elements(wave)

    openw, lun, outpath+'mike_thar_alpha_custom_shift.lst', /get_lun
    for ii = 0L, nall-1L do printf, lun, wave[ii], '1', 'ThAr', $
      format='(F9.4,1x,A1,1x,A4)'
    free_lun, lun
    
    
return
end
    
