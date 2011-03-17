pro extranucs_log12oh, hii, sings, result, write=write
; jm08feb20nyu - estimate metallicities for the extra-nuclear regions 

    path = sings_path(/projects)+'log12oh/'
    outpath = sings_path(/projects)+'extranucs/'
    version = sings_log12oh_version()

    nmonte = 500L
    
    if (n_elements(hii) eq 0L) then hii = $
      mrdfits(path+'sings_hiiregions_'+version+'.fits.gz',1,/silent)
    if (n_elements(sings) eq 0L) then sings = $
      mrdfits(path+'sings_log12oh_'+version+'.fits.gz',1,/silent)

    nucs = rsex('extranucs.sex')
    nnucs = n_elements(nucs)
    
    galaxynucs = strtrim(nucs[uniq(nucs.galaxy,sort(nucs.galaxy))].galaxy,2)

    result = replicate({radius: -999.0, rr25: -999.0, $
      log12oh_pt05: -999.0, log12oh_pt05_err: -999.0, $
      log12oh_kk04: -999.0, log12oh_kk04_err: -999.0},nnucs)
    result = struct_addtags(nucs,result)

    for ii = 0L, n_elements(galaxynucs)-1L do begin

       thesenucs = where((galaxynucs[ii] eq nucs.galaxy),nthese)
       match = where((galaxynucs[ii] eq strtrim(sings.galaxy,2)),nmatch)

       if sings[match].gradient_flag then begin
          
; deproject    
          
          nucra = 3600.0*15.0*hms2dec(sings[match].ra)
          nucdec = 3600.0*hms2dec(sings[match].dec)

          raoffset = (3600.0*15.0*im_hms2dec(nucs[thesenucs].ra)-nucra)*$
            cos(hms2dec(sings[match].dec)*!dtor)
          deoffset = 3600.0*im_hms2dec(nucs[thesenucs].dec)-nucdec

          result[thesenucs].radius = im_hiiregion_deproject(sings[match].incl,$ ; [arcsec]
            sings[match].pa,raoffset,deoffset)
          result[thesenucs].rr25 = result[thesenucs].radius/(sings[match].r25*60.0)

          coeff_pt05 = [sings[match].hii_pt05_log12oh_central[0],sings[match].hii_pt05_slope[0]]
          coeff_kk04 = [sings[match].hii_kk04_log12oh_central[0],sings[match].hii_kk04_slope[0]]
          coeff_pt05_err = [sings[match].hii_pt05_log12oh_central[1],sings[match].hii_pt05_slope[1]]
          coeff_kk04_err = [sings[match].hii_kk04_log12oh_central[1],sings[match].hii_kk04_slope[1]]
          
          result[thesenucs].log12oh_pt05 = poly(result[thesenucs].rr25,coeff_pt05)
          result[thesenucs].log12oh_kk04 = poly(result[thesenucs].rr25,coeff_kk04)

; Monte Carlo the errors          
          
          coeff_pt05_monte = fltarr(2,nmonte)
          coeff_pt05_monte[0,*] = coeff_pt05[0] + randomn(seed,nmonte)*coeff_pt05_err[0]
          coeff_pt05_monte[1,*] = coeff_pt05[1] + randomn(seed,nmonte)*coeff_pt05_err[1]
          log12oh_pt05_monte = fltarr(nmonte)
          for it = 0L, nthese-1L do begin
             for im = 0L, nmonte-1L do log12oh_pt05_monte[im] = poly(result[thesenucs[it]].rr25,coeff_pt05_monte[*,im])
             result[thesenucs[it]].log12oh_pt05_err = stddev(log12oh_pt05_monte)
          endfor
          
          coeff_kk04_monte = fltarr(2,nmonte)
          coeff_kk04_monte[0,*] = coeff_kk04[0] + randomn(seed,nmonte)*coeff_kk04_err[0]
          coeff_kk04_monte[1,*] = coeff_kk04[1] + randomn(seed,nmonte)*coeff_kk04_err[1]
          log12oh_kk04_monte = fltarr(nmonte)
          for it = 0L, nthese-1L do begin
             for im = 0L, nmonte-1L do log12oh_kk04_monte[im] = poly(result[thesenucs[it]].rr25,coeff_kk04_monte[*,im])
             result[thesenucs[it]].log12oh_kk04_err = stddev(log12oh_kk04_monte)
          endfor
          
       endif else begin

; not sure what to do here          
          
       endelse
       
    endfor

    if keyword_set(write) then begin

       txtoutfile = 'extranucs_log12oh.dat'
       splog, 'Writing '+outpath+txtoutfile+'.'
       openw, lun, outpath+txtoutfile, /get_lun
       printf, lun, '## Nebular Abundances for the SINGS Extranuclear HII Regions'
       printf, lun, '## John Moustakas, john.moustakas@nyu.edu, '+im_today()
;      printf, lun, '## NEARLY FINAL VALUES'
       printf, lun, '## PROVISIONAL'
;      printf, lun, '## '+version
       printf, lun, '## Discarded or "no measurement" are indicated using -999.'
       printf, lun, '#  1 Galaxy'
       printf, lun, '#  2 Region'
       printf, lun, '#  3 RA'
       printf, lun, '#  4 DEC'
       printf, lun, '#  5 Radius [deprojected, arcsec]'
       printf, lun, '#  6 R/R25 [deprojected]'
       printf, lun, '#  7 Log12oh_PT05'
       printf, lun, '#  8 Log12oh_PT05_Err'
       printf, lun, '#  9 Log12oh_KK04'
       printf, lun, '# 10 Log12oh_KK04_Err'
     
       struct_print, result, lun=lun, /no_head

    endif

return
end
