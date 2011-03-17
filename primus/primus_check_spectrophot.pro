pro primus_check_spectrophot, debug=debug, ps=ps
; jm08may12nyu - check the PRIMUS spectrophotometry

;   jm08jul24nyu - use the individual mask path!
    path = primus_path()

stop    
    
    if keyword_set(ps) then begin
       debug = 0L
       postthick1 = 4.0
       postthick2 = 3.0
       postthick3 = 8.0
    endif else begin
       if keyword_set(debug) then im_window, 0, $
         xratio=0.7, yratio=0.6
       postthick1 = 2.0
       postthick2 = 2.0
       postthick3 = 2.0
    endelse

; read the DEEP2 catalog and some k-correct preliminaries
    
    alldeep = read_deep2(/kcorr)
    vname = 'default.nolines'
    k_load_vmatrix, vmatrix, lambda, vname=vname

    bandindx = 1L ; normalize everything to the R-band

; some parameters of the comparison; only fit galaxies brighter than
; R=23.5; note that spectra of objects brighter than R~22 should have
; S/N>10

    zmin = 0.01
    zmax = 1.0    
    minmag = 22.0 ; 23.5 ; R-band

    plotscale = 1E17
    light = 2.99792458D18       ; speed of light [A/s]

; specify the rerun and mask(s)    
    
    rerun = '0213' ; "true" rerun = 0013

    allmasks = yanny_readone(getenv('PRIMUS_DIR')+'/data/primus_observed_masks.par')
    allmasks_det = yanny_readone(getenv('PRIMUS_DIR')+'/data/primus_observed_mask_details.par')
    these = where(strmatch(allmasks_det.field,'*deep*',/fold))
    mask_rootname = strtrim(allmasks_det[these].mask,2)
    mask_field = strtrim(allmasks_det[these].field,2)
    mask_date = strtrim(allmasks_det[these].date,2)
    nmask = n_elements(mask_rootname)
;   mask_rootname = 'd23a0045'
;   mask_rootname = 'vvds0138'

    good = lindgen(nmask)
    for imask = 0L, nmask-1L do begin
       thismask = getenv('PRIMUS_REDUX')+'/'+rerun+'/'+$
         mask_date[imask]+'/'+mask_rootname[imask]+'.extract.fits.gz'
;      print, imask, thismask
       if (file_test(thismask,/regular) eq 0L) then begin
          splog, 'Mask '+thismask+' not reduced...skipping.'
          good[imask] = -1L
       endif
    endfor

    good = good[where(good ne -1L,nmask)]
    mask_rootname = mask_rootname[good]
    mask_field = mask_field[good]
    mask_date = mask_date[good]

; initialize the output data structure

    result_template = {mask: '', field: '', date: '', $
      objno: 0L, z: -1.0, mag: 0.0, coeffs: fltarr(9), chi2: -1.0, $
      modelnorm: 0.0, primusnorm1: 0.0, primusnorm2: 0.0, $
      model_maggies_synth: fltarr(3), primus_maggies_synth1: fltarr(3), $
      primus_maggies_synth2: fltarr(3), maggies: fltarr(3)}
    
; use the defaults from PRIMUS_FIT_REDSHIFT 

    if (n_elements(ftweak) eq 0L) then ftweak = 0.05
    if (n_elements(order) eq 0L) then order = 3L
    if (n_elements(wavemin) eq 0L) then wavemin = 4500.
    if (n_elements(wavemax) eq 0L) then wavemax = 9500.

    constpsf = 1
    oiitweak = 0
    reweight = 0
    nointrachip = 0
    chooseslit = 1

; loop on each mask     
    
;   for imask = 15, nmask-1L do begin
    for imask = 0L, nmask-1L do begin

       qarootname = 'qa_spectrophot_'+mask_rootname[imask]+'_'+mask_date[imask]       
       
; read the 2D data and flag crappy spectra, again following
; PRIMUS_FIT_REDSHIFT; be careful with the NIGHT variable in this
; loop! 

       primus_read_1dinputs, mask_rootname[imask], rerun, ext=extract, $
         slits=slits, airmass=airmass, /nooned; night=mask_date[imask]
       photoinfo = primus_get_photomaggies(mask_rootname[imask],slits)

; reweight (whatever that means)
       if keyword_set(reweight) then begin
          for ii = 0, n_elements(extract) -1 do begin
             weight = construct_weight(extract[ii].wave1, ivar=extract[ii].fivar1, $
               newivar=newivar, flux=extract[ii].fopt1*0)
             extract[ii].fivar1 = newivar
             weight = construct_weight(extract[ii].wave2, ivar=extract[ii].fivar2, $
               newivar=newivar, flux=extract[ii].fopt2*0)
             extract[ii].fivar2 = newivar
          endfor
       endif
       
; following PRIMUS_FIT_REDSHIFT, do not fit to masked pixels
       extract.fivar1 = extract.fivar1 * (extract.mask1 eq 0)
       extract.fivar2 = extract.fivar2 * (extract.mask2 eq 0)
       pair = transpose(intarr(n_elements(extract), 2))
       pair[0,*] = 1
       pair[1,*] = 2

       if keyword_set(chooseslit) then begin
          for i = 0L, n_elements(extract)-1L do begin
             IF (extract[i].badextract AND primus_flagval('SPEXTRACT', 'A_FULLMASK')) GT 0 THEN $
               pair[0,i] = 0
             IF (extract[i].badextract AND primus_flagval('SPEXTRACT', 'B_FULLMASK')) GT 0 THEN $
               pair[1,i] = 0
          endfor
       endif

; only fit objects with existing DEEP2 spectroscopic redshifts and for
; which both slit A&B spectra are good (this code is lifted from
; PRIMUS_FIT_REDSHIFT, corresponding to ZCAL=1); note:

;    CURRZ - best estimate of redshift from external source
;    ZTYPE - provenance of the redshift:
;                   0 - SDSS
;                   1,6 - DEEP2
;                   2 - VVDS
;                   3 - COMBO-17 photo-z
;                   4 - USNO (assumes it is a star)
;                   5 - is an SDSS star
;                   8 - QSO photo-z

       index = where(((slits.ztype eq 1) or (slits.ztype eq 6)) and $
         (slits.currz gt zmin) and (slits.currz lt zmax) and $
         (pair[0,*] gt 0) and (pair[1,*] gt 0) and (slits.mag lt minmag))
       keep = where(total((extract[index].wave1 gt 6000) and (extract[index].wave1 lt 10000), 1) gt 10 and $
         total((extract[index].wave2 gt 6000) and (extract[index].wave2 lt 10000), 1) gt 10)
       index = index[keep] & nindex = n_elements(index)

; cross-match against the parent DEEP2 catalog and store some info;
; use a FOR loop to preserve the order

;      spherematch, alldeep.ra, alldeep.dec, slits[index].ra, $
;        slits[index].dec, 3.0/3600.0, m1, m2
       deepmatch = lonarr(nindex)
       for ii = 0L, nindex-1L do deepmatch[ii] = $
         where(slits[index[ii]].objno eq alldeep.objname)
       flag = where(deepmatch eq -1L,nflag,comp=good)
       if (nflag ne 0L) then begin
;         message, 'This should not happen!'
          for iflag = 0L, nflag-1L do $
            splog, 'Object(s) '+string(slits[index[flag[iflag]]].objno,$
            format='(I0)')+' not found in my DEEP2 catalog!'
          deepmatch = deepmatch[good]
          index = index[good]
          nindex = n_elements(good)
       endif
       deep = alldeep[deepmatch]

       result = replicate(result_template,nindex)
       result.mask = mask_rootname[imask]
       result.field = mask_field[imask]
       result.date = mask_date[imask]
       result.objno = deep.objname

; open up the PS file and then loop on each object       
       
       if keyword_set(ps) then begin
          splog, 'Writing '+path+qarootname+'.ps'
          dfpsplot, path+qarootname+'.ps', /landscape, /color
       endif

       for iobj = 4L, nindex-1L do begin
;      for iobj = 0L, nindex-1L do begin

          print, format='("Galaxy ",I0,"/",I0,A10,$)', $
            iobj+1, nindex, string(13b)
          zz = slits[index[iobj]].currz
          filters = photoinfo[index[iobj]].filterlist
          filtinfo = im_filterspecs(filterlist=filters)

; fit the PRIMUS spectra at the known redshift          
          
          chi2 = primus_k_fit_spec(extract[index[iobj]],zz,airmass=airmass,$
            pair=pair[*,index[iobj]],constpsf=constpsf,grid=grid, $
            rcounts=rcounts,tweak=tweak,wavemin=wavemin,wavemax=wavemax,$
            order=order,ftweak=ftweak,nofflux=nofflux,nointrachip=nointrachip,$
            zwarning=gzwarning,primus=primus,sky=sky,oiitweak=oiitweak,$
            coeffs=coeffs,highwave=bestwave,highflux=highflux,$
            photoinfo=photoinfo[index[iobj]],photoz=1)
          if (finite(chi2) eq 0) then begin
             splog, 'Infinite chi2 on object (index) ', index[iobj]
             result[iobj].chi2 = -999.0
             continue
          endif

; re-build the best-fitting high-resolution template; note that the
; factor of (1+z) has not been applied to BESTWAVE, but it *has* to
; HIGHFLUX (see PRIMUS_GALAXY_TEMPLATES); there is some code here to
; crop to [4000,9500], but don't do it because we want to be
; able to extrapolate 

          bestwave = bestwave*(1.0+zz)
          bestflux = highflux#coeffs
;         use_bestwave = where((bestwave ge wavemin) and (bestwave le wavemax))
;         bestwave = bestwave[use_bestwave]
;         bestflux = bestflux[use_bestwave]
          bestwave_edges = k_lambda_to_edges(bestwave)

;         bestflux_flam = highflux#coeffs
;         bestflux_fnu1 = bestflux_flam*bestwave^2.0/light
;         bestflux_fnu2 = (highflux*rebin(reform(bestwave,n_elements(bestwave),1),$
;           n_elements(bestwave),9)^2.0/light)#coeffs
;         djs_plot, bestwave, bestflux_fnu1, ps=10, xsty=3, ysty=3, charsize=2
;         djs_oplot, bestwave, bestflux_fnu2, ps=10, color='red'
          
; convert the primus spectra to physical units (erg/s/cm2/A)

          wave1 = extract[index[iobj]].wave1
          wave2 = extract[index[iobj]].wave2
          counts2flam1 = template2primus(wave1,wave1*0.0+1.0,$
            airmass=airmass,slitinfo=extract[index[iobj]],pair=1)
          counts2flam2 = template2primus(wave2,wave2*0.0+1.0,$
            airmass=airmass,slitinfo=extract[index[iobj]],pair=2)
          use_wave1 = where((wave1 ge wavemin) and (wave1 le wavemax))
          use_wave2 = where((wave2 ge wavemin) and (wave2 le wavemax))
          wave1 = wave1[use_wave1] & wave2 = wave2[use_wave2]

          flux1 = extract[index[iobj]].fopt1[use_wave1]/counts2flam1[use_wave1]
          flux2 = extract[index[iobj]].fopt2[use_wave2]/counts2flam2[use_wave2]
          model_flux1 = rcounts[use_wave1,0]/counts2flam1[use_wave1]
          model_flux2 = rcounts[use_wave2,1]/counts2flam2[use_wave2]
          
;         counts2flam = primus_fluxvector(extract[index[iobj]])
;         flux1 = extract[index[iobj]].fopt1[use_wave1]/counts2flam[use_wave1,0]*$
;           extract[index[iobj]].calib1[use_wave1]
;         flux2 = extract[index[iobj]].fopt2[use_wave2]/counts2flam[use_wave2,1]*$
;           extract[index[iobj]].calib2[use_wave2]
;         model_flux1 = rcounts[use_wave1,0]/counts2flam[use_wave1,0]*$
;           extract[index[iobj]].calib1[use_wave1]
;         model_flux2 = rcounts[use_wave2,1]/counts2flam[use_wave2,1]*$
;           extract[index[iobj]].calib2[use_wave2]
          
;         thru = primus_throughput(extract[index[iobj]].wave1,$
;           ccd=extract[index[iobj]].ccdnum,airmass=airmass)

; synthesize broadband photometry based on the best-fitting template;
; also synthesize photometry directly from the spectra, but note that
; there may be extrapolation issues
          
          model_maggies_synth = reform(k_project_filters(bestwave_edges,bestflux, $
            filterlist=filters,/silent))
          primus_maggies_synth1 = reform(k_project_filters(k_lambda_to_edges(wave1),$
            flux1,filterlist=filters,/silent)) ; could use MODEL_FLUX1
          primus_maggies_synth2 = reform(k_project_filters(k_lambda_to_edges(wave2),$
            flux2,filterlist=filters,/silent)) ; could use MODEL_FLUX2

; scale all the spectra relative to the R-band
          
          modelnorm = photoinfo[index[iobj]].maggies[bandindx]/model_maggies_synth[bandindx]
          primusnorm1 = photoinfo[index[iobj]].maggies[bandindx]/primus_maggies_synth1[bandindx]
          primusnorm2 = photoinfo[index[iobj]].maggies[bandindx]/primus_maggies_synth2[bandindx]

; store the results          
          
          result[iobj].z = zz
          result[iobj].chi2 = chi2
          result[iobj].mag = slits[index[iobj]].mag ; R-band, for DEEP2
          result[iobj].modelnorm = modelnorm
          result[iobj].primusnorm1 = primusnorm1
          result[iobj].primusnorm2 = primusnorm2
          result[iobj].coeffs = coeffs
          result[iobj].model_maggies_synth = model_maggies_synth*modelnorm
          result[iobj].primus_maggies_synth1 = primus_maggies_synth1*primusnorm1
          result[iobj].primus_maggies_synth2 = primus_maggies_synth2*primusnorm2
          result[iobj].maggies = photoinfo[index[iobj]].maggies
;         print, maggies_synth/photoinfo[index[iobj]].maggies

; debugging plot          

          if keyword_set(debug) or keyword_set(ps) then begin

; first plot the k-correct fit to the broadband photometry             
             
             fnu2flam = 10^(-0.4*48.6)*light/photoinfo[index[iobj]].weff^2.0
             wave = k_lambda_to_centers(lambda)*(1.0+zz)
             flux_flam = plotscale*(vmatrix#deep[iobj].coeffs)/(1.0+zz)
             maggies_flam = plotscale*photoinfo[index[iobj]].maggies*fnu2flam
             maggies_flam_ivar = photoinfo[index[iobj]].maggieivar/fnu2flam^2.0/(plotscale)^2.0

             xrange = [3000,10000] & get_element, wave, xrange, xx
             yrange = fltarr(2)
             yrange[0] = 0.0 ; NOTE!
;            yrange[0] = min(flux_flam[xx[0]:xx[1]])<min(plotscale*bestflux*modelnorm)
;            yrange[1] = max(flux_flam[xx[0]:xx[1]])>max(plotscale*bestflux*modelnorm)
;            yrange[0] = min(flux_flam[xx[0]:xx[1]])<min(plotscale*bestflux*modelnorm)<$
;              min(flux1*primusnorm)<min(flux2*primusnorm)
             ss1 = im_stats(plotscale*model_flux1*primusnorm1,sigrej=3.0)
             ss2 = im_stats(plotscale*model_flux2*primusnorm2,sigrej=3.0)
             ss3 = im_stats(flux_flam[xx[0]:xx[1]],sigrej=3.0)
             ss4 = im_stats(plotscale*bestflux*modelnorm,sigrej=3.0)
             yrange[1] = ss1.maxrej>ss2.maxrej>ss3.maxrej>ss4.maxrej
;            yrange[1] = max(plotscale*bestflux*modelnorm)

             djs_plot, [0], [0], /nodata, charsize=1.8, xthick=postthick1, ythick=postthick1, $
               title='DEEP2/'+string(result[iobj].objno,format='(I0)')+$
               ', z = '+strtrim(string(result[iobj].z,format='(F12.4)'),2)+$
               ', R = '+strtrim(string(result[iobj].mag,format='(F12.2)'),2), $
               xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
               xtitle='Observed Wavelength (\AA)', $
               ytitle='Flux (10^{-17} '+flam_units()+')', charthick=postthick2
             djs_oplot, wave, flux_flam, ps=10, color='grey', thick=postthick1
             plotsym, 8, 3.0, fill=1, thick=postthick3
             djs_oplot, photoinfo[index[iobj]].weff, plotscale*fnu2flam*result[iobj].maggies, $
               ps=8, color=djs_icolor('dark green'), thick=postthick1

; now plot the best-fitting (high-resolution) spectral template, the
; photometry synthesized from the best-fit template, and the actual
; PRIMUS spectra

;            put this plot at the beginning of the individual spectral
;            plot!
stop             
             
             djs_oplot, bestwave, plotscale*bestflux*modelnorm, color='blue', thick=postthick1
             plotsym, 0, 3.0, fill=0, thick=postthick3
             djs_oplot, photoinfo[index[iobj]].weff, plotscale*fnu2flam*result[iobj].model_maggies_synth, $
               ps=8, color=djs_icolor('orange'), thick=postthick1

             djs_oplot, wave1, plotscale*flux1*primusnorm1, line=1, $
               color='purple', thick=postthick1
             djs_oplot, wave2, plotscale*flux2*primusnorm2, line=1, $
               color='red', thick=postthick1
             djs_oplot, wave1, plotscale*model_flux1*primusnorm1, line=0, thick=postthick2
             djs_oplot, wave2, plotscale*model_flux2*primusnorm2, line=0, thick=postthick2

             for ifilt = 0L, n_elements(filters)-1L do begin
                fthese = lindgen(filtinfo[ifilt].filtn)
                djs_oplot, filtinfo[ifilt].filtw[fthese], filtinfo[ifilt].filtf[fthese]/$
                  max(filtinfo[ifilt].filtf[fthese])*0.2*(yrange[1]-yrange[0]), $
                  line=1, thick=postthick2
             endfor
             
             label = ['K-correct fit to BRI','PRIMUS best-fit model',$
               'Slit A','Slit B']
;            label = ['K-correct fit to BRI','PRIMUS best-fit model '+$
;              '(\chi^{2} = '+strtrim(string(chi2,format='(F12.1)'),2)+')',$
;              'Slit A','Slit B']
             legend, textoidl(label), /left, /top, box=0, clear=keyword_set(ps), $
               charsize=1.6, charthick=postthick2, thick=postthick3, line=[0,0,1,1], $
               color=djs_icolor(['grey','blue','purple','red'])
;            legend, ['BRI photometry','k-correct best-fit','PRIMUS best fit',$
;              'Slit A','Slit B'], /left, /top, box=0, /clear, charsize=1.6, $
;              charthick=postthick2, thick=postthick1, line=[0,0,0,1,1], $
;              color=djs_icolor(['dark green','grey','blue','orange','red'])

;            label = ['BRI - DEEP2','BRI - Synthesized']
;            thisthick = postthick3
;            im_legend, textoidl(label), /right, /bottom, box=0, /clear, charsize=1.6, $
;              charthick=postthick2, fill=[1,0], thick=thisthick, $
;              symsize=1.5, spacing=1.5, $
;              color=djs_icolor(['dark green','orange']), psym=[108,106]
             
             if keyword_set(debug) then cc = get_kbrd(1)

          endif
             
       endfor 

       if keyword_set(ps) then begin
          dfpsclose
          spawn, 'gzip -f '+path+qarootname+'.ps'
       endif
       splog, 'Writing '+path+qarootname+'.fits'
       mwrfits, result, path+qarootname+'.fits', /create

    endfor ; close the MASK loop
          
stop    
    
return
end
