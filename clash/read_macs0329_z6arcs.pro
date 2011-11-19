function addphot, in, all
    filt = clash_filterlist(short=short,zpt=zpt)
    out = in
    nall = n_elements(all)
    for ii = 0, n_elements(short)-1 do begin
       gd = where(finite(all.(tag_indx(all,short[ii]+'_fluxerr'))))
       totflux = total(all[gd].(tag_indx(all[gd],short[ii]+'_flux')))
       totferr = sqrt(total(all[gd].(tag_indx(all[gd],short[ii]+'_fluxerr'))^2))
       if (totflux le 0) then begin
          out.(tag_indx(out,short[ii]+'_mag')) = 99.0 
          out.(tag_indx(out,short[ii]+'_magerr')) = -2.5*alog10(totferr)+zpt[ii]
       endif else begin
          out.(tag_indx(out,short[ii]+'_mag')) = -2.5*alog10(totflux)+zpt[ii]
          out.(tag_indx(out,short[ii]+'_magerr')) = 2.5*totferr/totflux/alog(10)
       endelse
    endfor
return, out
end

function read_macs0329_z6arcs, adi=adi, usemag=usemag
; jm11nov11ucsd - read the catalog photometry

CODE RELEGATED BY READ_Z6ARCS_PHOTOMETRY!
    
    isedpath = clash_path(/ised)
    datapath = clash_path(/macs0329_z6arcs)

    filt = clash_filterlist(short=short,weff=weff,zpt=zpt,nice=nice)
    nfilt = n_elements(filt)
    
    adi = rsex(datapath+'macs0329_z6arcs.cat')
    struct_print, adi
    ngal = n_elements(adi)

    irac = rsex(datapath+'macs0329_z6arcs_irac.cat')
;   irac = rsex(datapath+'macs0329_z6arcs_irac_larry.cat')
;   irac = rsex(datapath+'macs0329_z6arcs_irac_leonidas.cat')

    if keyword_set(usemag) then begin

       allcat = read_clash_catalog('macs0329',/ir)

;; old code:    
;       spherematch, allcat.ra, allcat.dec, 15D*hms2dec(adi[3].ra), $
;         hms2dec(adi[3].dec), 1D/3600, m1, m2, maxmatch=0
;       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
;       cat = allcat[m1]
       
       cat = {id: 0L, ra: 0.0D, dec: 0.0D, ch1_mag: -99.0, ch1_magerr: -99.0, $
         ch2_mag: -99.0, ch2_magerr: -99.0}
       for ii = 0, nfilt-1 do cat = create_struct(cat,$
         short[ii]+'_mag',-99.0,short[ii]+'_magerr',-99.0)
       cat = replicate(cat,ngal)
       
; quick work-around for shredded arcs; just add photometry linearly
; and the errors in quadrature

; arc 1.1 = 988, 989, 990; CTE problem on F336W, replace with upper
; limit  
       these = where(allcat.id ge 988 and allcat.id le 990)
       cat[0] = addphot(cat[0],allcat[these])
       cat.f336w_magerr = cat.f336w_mag
       cat.f336w_mag = 99.0
       
; arc 1.2 = 1132, 1133, 1134, 1135
       these = where(allcat.id ge 1132 and allcat.id le 1135)
       cat[1] = addphot(cat[1],allcat[these])
       
; arc 1.2 = 1388, 1387
       these = where(allcat.id ge 1387 and allcat.id le 1388)
       cat[2] = addphot(cat[2],allcat[these])
       
; arc 1.4 = 517 - ok?
       these = where(allcat.id eq 517)
       cat[3] = addphot(cat[3],allcat[these])

; add IRAC
       minerr = 0.0
       cat.ch1_mag = irac.ch1_ab
       cat.ch2_mag = irac.ch2_ab
       cat.ch1_magerr = sqrt(irac.ch1_aberr^2+minerr^2)
       cat.ch2_magerr = sqrt(irac.ch2_aberr^2+minerr^2)

    endif else begin
; use Marc's manual photometry of the arc
       cat = {id: 0L, ra: 0.0D, dec: 0.0D, ch1_flux: -99D, ch1_fluxerr: -99D, $
         ch2_flux: -99D, ch2_fluxerr: -99D}
       for ii = 0, nfilt-1 do cat = create_struct(cat,$
         short[ii]+'_flux',-99D,short[ii]+'_fluxerr',-99D)
       cat = replicate(cat,ngal)

       file = datapath+'z6arc_photometry_obj1_'+['1','2','3','4']+'.out'
       for jj = 0, ngal-1 do begin
          readcol, file[jj], obj, flux, fluxerr, mag, magerr, hstzpt, $
            hstzptext, wave, filename, format='X,A,F,F,F,F,F,F,F,A', $
            comment='#', /silent
          get_element, weff, wave, windx
          for ii = 0, n_elements(windx)-1 do begin
             fluxindx = tag_indx(cat,short[windx[ii]]+'_flux')
             fluxerrindx = tag_indx(cat,short[windx[ii]]+'_fluxerr')
             cat[jj].(fluxindx) = flux[ii]
             cat[jj].(fluxerrindx) = fluxerr[ii]

; totally ad hoc!!!!
             if strmatch(nice[windx[ii]],'*uvis*',/fold) then cat[jj].(fluxindx) = 0.0
;            if weff[windx[ii]]/7.18 le 1215.0 then cat[jj].(fluxindx) = 0.0
;            print, short[windx[ii]], weff[windx[ii]], wave[ii], mag[ii], $
;              -2.5*alog10(flux[ii])+hstzpt[ii], -2.5*alog10(flux[ii])+zpt[windx[ii]]
          endfor
       endfor

; add IRAC
       minerr = 0.0
       cat.ch1_flux = 10D^(-0.4*irac.ch1_ab)
       cat.ch2_flux = 10D^(-0.4*irac.ch2_ab)
       cat.ch1_fluxerr = (alog(10)/2.5)*sqrt(irac.ch1_aberr^2+minerr^2)*cat.ch1_flux
       cat.ch2_fluxerr = (alog(10)/2.5)*sqrt(irac.ch2_aberr^2+minerr^2)*cat.ch2_flux

;      lim = where(irac.ch2 gt 90.0)
;      cat[lim].ch2_flux = 0D
;      cat[lim].ch2_fluxerr = 10D^(-0.4*irac[lim].ch2_err)/2.0 ; 1-sigma
    endelse
       
return, cat
end

;; photometry from Leonidas:
;    cat.ch1_mag = -2.5*alog10(irac.ch1*10.0^(-0.4*23.9))
;    cat.ch1_magerr = 2.5*irac.ch1_err/irac.ch1/alog(10)
;
;    gd = where(irac.ch2 gt 0.0)
;    cat[gd].ch2_mag = -2.5*alog10(irac[gd].ch2*10.0^(-0.4*23.9))
;    cat[gd].ch2_magerr = 2.5*irac[gd].ch2_err/irac[gd].ch2/alog(10)
;
;    lim = where(irac.ch2 lt 0.0 and irac.ch2_err gt 0.0)
;    cat[lim].ch2_mag = 99.0
;    cat[lim].ch2_magerr = -2.5*alog10(irac[lim].ch2_err*10.0^(-0.4*23.9))

