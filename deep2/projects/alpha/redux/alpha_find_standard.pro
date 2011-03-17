function alpha_find_standard, stdinfo, search=search
; jm09jan21nyu - return the appropriate ESOFIL; add more standards
; from here:
; http://www.ctio.noao.edu/spectrographs/4m_Echelle/standards.html

    if (n_elements(search) eq 0L) then search = 60.0 ; 1 arcmin

    nstd = n_elements(stdinfo)
    if (nstd gt 1L) then begin
       for jj = 0L, nstd-1L do begin
          info1 = alpha_find_standard(stdinfo[jj],search=search)
          if (jj eq 0L) then info = info1 else info = [info,info1]
       endfor
       return, info
    endif

    nstar = 3
    std = {hr: '', name: '', stdra: 0.0D, stddec: 0.0D, esofil: ''}
    std = replicate(std,nstar)

; HR3454    
    std[0].hr      = 'HR3454'
    std[0].name    = 'eta Hya'
    std[0].stdra   = 15.0D*hms2dec('08:43:13.4')
    std[0].stddec  = hms2dec('03:23:55')
    std[0].esofil  = 'fhr3454.dat'
; HR4468
    std[1].hr      = 'HR4468'
    std[1].name    = 'theta Crt'
    std[1].stdra   = 15.0D*hms2dec('11:36:40.8')
    std[1].stddec  = hms2dec('-09:48:08')
    std[1].esofil  = 'fhr4468.dat'
; HR7950
    std[2].hr      = 'HR7950'
    std[2].name    = 'epsilon Aqr'
    std[2].stdra   = 15.0D*hms2dec('20:47:40.5')
    std[2].stddec  = hms2dec('-09:29:45')
    std[2].esofil  = 'fhr7950.dat'

; these are not actual standard stars, but we're going to try
; and use HR7950     
    if strmatch(strcompress(stdinfo.obj,/remove),'*HR6930*',/fold) or $
      strmatch(strcompress(stdinfo.obj,/remove),'*HR7830*',/fold) then begin
       m1 = 2L 
    endif else begin
       spherematch, std.stdra, std.stddec, stdinfo.ra, stdinfo.dec, $
         search/3600.0, m1, m2
       if (m1[0] eq -1L) then message, 'Unrecognized standard star'
    endelse

return, std[m1]
end
