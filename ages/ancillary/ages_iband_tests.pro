pro ages_iband_tests
; jm10feb02ubsd - tests of the new and old I-band photometry

    common iband_tests, phot, diff, mag

    if (n_elements(phot) eq 0) then begin
       photfile = ages_path(/mycatalogs)+'ages_photometry_'+$
         ages_version(/photometry)+'.fits.gz'
       splog, 'Reading '+photfile
       phot = mrdfits(photfile,1)
       phot = phot[where(phot.z gt 0.0)]
    endif

    if (n_elements(mag) eq 0) then begin
       usnofile = ages_path(/mycatalogs)+'ages_usno.fits.gz'
       usno = mrdfits(usnofile,1)
       keep = where(usno.rmag lt 12.0,nkeep)
       usno = usno[keep]

       ing = spheregroup(usno.ra,usno.dec,30.0/3600.0,$
         firstg=firstg,multg=multg,nextg=nextg)
       firstg = firstg[0L:max(ing)]
       usno = usno[firstg]
       nstar = n_elements(usno)
       
       mag = fltarr(n_elements(phot))
       diff = fltarr(n_elements(phot))
       for ii = 0L, n_elements(phot)-1L do begin
          stardiff = djs_diff_angle(phot[ii].ra,phot[ii].dec,usno.ra,usno.dec)
          diff[ii] = min(stardiff,indx)
          mag[ii] = usno[indx].rmag
       endfor
    endif

; plot the difference in the I(Kron) and I(Kron,R) magnitudes versus
; the distance to the nearest USNO star
    
    bstar = where(phot.bstar,comp=clean)

    ircolor = phot.i_mag_aper_06-phot.r_mag_aper_06
    ikron = phot.i_auto            ; I_Kron
    irband = phot.r_auto + ircolor ; I_R = R_Kron+(I-R)_6"
    
    dmag = max(mag[clean])-mag[clean]
    symsize = 0.8*10.0^(dmag/max(dmag))
;   symsize = alog10(10.0*mag[clean]/min(mag[clean]))

    psfile = ages_path(/qaplots)+'ages_brightstar_tests.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.6
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
      xr=[0.1,10], yr=[-2.5,2.5], ytitle='I(Kron)-I(Kron,R) (mag)', $
      xtitle='Distance to Nearest USNO Star (arcmin)', /xlog
    djs_oplot, diff[clean]*60.0, ikron[clean]-irband[clean], $
      psym=6, symsize=symsize
;   djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
;     xr=[0.1,10], yr=[-1,1], ytitle='I(Kron)-I(Kron,R) (mag)', $
;     xtitle='Distance to Nearest USNO Star (arcmin)', /xlog
;   djs_oplot, diff[clean]*60.0, ikron[clean]-irband[clean], $
;     psym=6, symsize=symsize
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
stop    
    
; reproduce Fig 2 (the bright-star problem plot) in the 2005
; documentation
    bstar = where(phot.bstar,comp=clean)
    ikron = phot.i_mag_auto
    iaper = phot.i_mag_aper_06

    psfile = ages_path(/qaplots)+'ages_iband_tests.ps'
    im_plotconfig, 6, pos, psfile=psfile, charsize=1.6
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos[*,0], $
      xr=[10,20], yr=[-0.5,5], xtickname=replicate(' ',10), $
      ytitle='I(6" aper)-I(Kron) (mag)'
    djs_oplot, ikron[clean], iaper[clean]-ikron[clean], psym=6, sym=0.2
    legend, 'brightstar=0', /left, /bottom, box=0
    
    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, position=pos[*,1], $
      xr=[10,20], yr=[-0.5,5], ytitle='I(6" aper)-I(Kron) (mag)', $
      xtitle='I(Kron) (mag)'
    djs_oplot, ikron[bstar], iaper[bstar]-ikron[bstar], psym=6, sym=0.2
    legend, 'brightstar=1', /left, /bottom, box=0
    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    
    
; --------------------------------------------------
; compare the new and old I-band photometry
    diff = abs(phot.i_auto-phot.i_mag_auto)
    bad = where(diff gt 0.01,comp=good)
    niceprint, phot[bad].i_auto, phot[bad].i_mag_auto, phot[bad].i_flag_subfield
    out = im_struct_trimtags(phot[bad],select=['i_alpha_j2000','i_delta_j2000',$
      'i_auto','i_mag_auto','i_mag_aper_10','i_flag_subfield'],$
      newtags=['ra','dec','i_auto_old','i_auto_new','i_aper_10','flag_subfield'])
    struct_print, out
    outfile = '~/iband_tests.fits'
    im_mwrfits, out, outfile
    
stop

return
end
    
