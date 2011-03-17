pro write_calib_planfile, planfile, calibplanfile
; write out a plan file for just the calibration data
    
    plan = yanny_readone(planfile,hdr=hdr)
    calib = where(strtrim(plan.flavor,2) ne 'science')
    new = plan[calib]

    plotfile = where(strmatch(hdr,'*plotfile*'))
    hdr[plotfile] = 'plotfile '+repstr(calibplanfile,'.par','.ps')
    logfile = where(strmatch(hdr,'*logfile*'))
    hdr[logfile] = 'logfile '+repstr(calibplanfile,'.par','.log')

    splog, 'Writing '+calibplanfile
    yanny_write, calibplanfile, ptr_new(new), hdr=hdr, /align

return
end

pro atlas_edit_plan, night, planfile, apinfofile, sings=sings, atlas=atlas
; each plan file must be custom-edited
    case night of
       '01nov10': newplan = newplan_01nov10(planfile,hdr=hdr,sings=sings,atlas=atlas,apinfo=apinfo)
       '01nov11': newplan = newplan_01nov11(planfile,hdr=hdr,sings=sings,atlas=atlas,apinfo=apinfo)
       '01nov13': newplan = newplan_01nov13(planfile,hdr=hdr,sings=sings,atlas=atlas,apinfo=apinfo)
       else: stop
    endcase
; write out

    splog, 'Writing '+planfile
    yanny_write, planfile, ptr_new(newplan), hdr=hdr, /align

    splog, 'Writing '+apinfofile
    im_mwrfits, apinfo, apinfofile, /clobber
    
return
end
    

function atlas_get_nights, run
; given the run, return the list of nights to reduce
    case run of
       '94nov': allnight = ''
       '95mar': allnight = ''
       '95oct': allnight = ''
       '96apr': allnight = ''
       '97apr': allnight = ''
       '98apr': allnight = ''
       '98jun': allnight = ''
       '98mar': allnight = ''
       '98oct': allnight = ''
       '99apr': allnight = ''
       '99may': allnight = ''
       '99nov': allnight = ''
       '00apr': allnight = ''
       '01nov': allnight = ['01nov10']
;      '01nov': allnight = ['01nov10','01nov11','01nov13']
       '01dec': allnight = ''
       '02apr': allnight = ''
       '02feb': allnight = ''
       '02may': allnight = ''
       '03may': allnight = ''
       '05apr': allnight = ''
       '06mar': allnight = ''
       '06may': allnight = ''
       else: message, 'Run not recognized'
    endcase
return, allnight
end

pro atlas_reduce_all, datapath, date, atlas=atlas, sings=sings, $
  plan=plan, calib=calib, standards=standards, sensfunc=sensfunc, $
  science=science, ccdproc=ccdproc, extract=extract, coadd=coadd, $
  fluxcal=fluxcal, final_spec1d=final_spec1d, clobber=clobber
; jm09dec17ucsd - uber-wrapper to reduce *all* the atlas/sings data 

; atlas_reduce_all, /plan, /calib, /stand, /clobber    
    
    atlaspath = getenv('RESEARCHPATH')+'/data/atlas/'
    spec1dpath = atlaspath+'spec1d/'

    sings = 1 & atlas = 0
;   sings = 0 & atlas = 1
    
    if keyword_set(sings) then begin
       splog, 'Reducing SINGS data'
;      allruns = ['01nov','01dec','02feb','02may','05apr','06mar','06may']
       allruns = ['01nov']
       suffix = 'sings'
    endif else begin
       splog, 'Reducing ATLAS data'
;      allruns = ['94nov','95mar','95oct','96apr',$
;        '97apr','98apr','98jun','98mar','98oct','99apr','99may',$
;        '99nov','00apr','01nov','02apr','02feb','02may','03may']
       allruns = ['01nov']
       suffix = 'atlas'
    endelse
    nrun = n_elements(allruns)

    for irun = 0, nrun-1 do begin
       rundir = atlaspath+allruns[irun]+'/'
       badpixfile = rundir+'badpix.dat'
       sensfuncfile = rundir+'sensfunc_'+allruns[irun]+'.fits'

       allnight = atlas_get_nights(allruns[irun])
       nnight = n_elements(allnight)

       for inight = 0, nnight-1 do begin
          datapath = rundir+allnight[inight]+'/'
          calibplanfile = datapath+'plan_calib_'+allnight[inight]+'.par'
          planfile = datapath+'plan_'+suffix+'_'+allnight[inight]+'.par'
          apinfofile = datapath+'apinfo_'+suffix+'_'+allnight[inight]+'.fits'
; --------------------------------------------------
; make the plan files
          if keyword_set(plan) then begin
             long_plan, '*.fits.gz', datapath+'Raw/', $
               planfile=planfile
             atlas_edit_plan, allnight[inight], planfile, $
               apinfofile, sings=sings, atlas=atlas
          endif

; --------------------------------------------------
; reduce the calibration data
          if keyword_set(calib) then begin
             if keyword_set(clobber) then write_calib_planfile, $
               planfile, calibplanfile
             pushd, datapath
             atlas_long_reduce, calibplanfile, /justcalib, $
               /verbose, calibclobber=clobber
             popd
          endif

; --------------------------------------------------
; reduce and extract the standards
          if keyword_set(standards) then begin
             pushd, datapath
             atlas_long_reduce, calibplanfile, badpixfile=badpixfile, $
               /juststd, /novac, clobber=clobber
             popd
          endif

; --------------------------------------------------
; reduce all the galaxies
          if keyword_set(science) then begin
             pushd, datapath
             atlas_long_reduce, planfile, badpixfile=badpixfile, $
               /justsci, /novac, skytrace=0, clobber=clobber, /chk
             popd
          endif

;; --------------------------------------------------
;; apply the basic calibration steps to the science data
;          if keyword_set(ccdproc) then begin
;             pushd, datapath
;             atlas_long_reduce, planfile, /justsci, /ccdproc, $
;               badpixfile=badpixfile, clobber=clobber
;             popd
;          endif

; --------------------------------------------------
; stitch E/W spectra together!
                    
          
; --------------------------------------------------
; extract the individual spectra
          if keyword_set(extract) then begin
             pushd, datapath
             atlas_long_extract, planfile, /justsci, $
               /ccdproc, clobber=clobber
             popd
          endif

      endfor ; NIGHT

; --------------------------------------------------
; make the sensitivity function using the standards from all the
; nights 
       if keyword_set(sensfunc) then begin
          for inight = 0, nnight-1 do begin
             datapath = atlaspath+allruns[irun]+'/'+allnight[inight]+'/'
             planfile = datapath+'plan_'+allnight[inight]+'.par'
             plan1 = yanny_readone(planfile)
             plan1.filename = datapath+'Science/std-'+strtrim(plan1.filename,2)
             if (inight eq 0) then plan = plan1 else $
               plan = [plan,plan1]
          endfor
          std = where((strmatch(plan.starname,'*...*') eq 0),nstd)
          splog, 'Found '+string(nstd,format='(I0)')+' standards'
          stdfiles = strtrim(plan[std].filename,2)
          std_names = strtrim(plan[std].starfile,2)
          
          sens = atlas_long_sensfunc(stdfiles,sensfuncfile,$
            std_name=std_names,/msk_balm)
       endif
       
; --------------------------------------------------
; coadd sequential observations of the same object; to select the list
; of files we read in the ispec crsplit files
       if keyword_set(coadd) then begin
          for inight = 0, nnight-1 do begin
             datapath = atlaspath+allruns[irun]+'/'+allnight[inight]+'/'
; coadd the CRsplits
             coaddfile = datapath+'crsplits_'+allnight[inight]+'.txt'
             if file_test(coaddfile,/reg) then begin
                coaddtxt = djs_readlines(coaddfile)
                coaddtxt = coaddtxt[where((strcompress(coaddtxt,/remove) ne '') and $
                  strmatch(coaddtxt,'*#*') eq 0)]
                message, 'if the file does not exist then doom!', /inf
                for ii = 0, n_elements(coaddtxt)-1 do begin
                   infiles = datapath+'Science/'+strtrim($
                     strsplit(coaddtxt[ii],' ',/extract),2)
                   infiles = repstr(infiles,'ra.','sci-a.')+'.gz' ; convert ispec file name
                   nobs = n_elements(infiles)
                   outfile = atlaspath+allruns[irun]+'/spec1d_unfluxed/'+$
                     file_basename(strmid(infiles[0],0,strpos(infiles[0],'.fits.gz'))+$
                     '_'+string(nobs,format='(I0)')+'.fits')
                   outfile = repstr(outfile,'sci-','scr') ; hack!
                   long_coadd, infiles, 1, outfil=outfile
                endfor 
             endif
; now deal with single exposures
             stop
          endfor 
       endif 

; --------------------------------------------------
; flux-calibrate
       if keyword_set(fluxcal) then begin
          allfiles = file_search(atlaspath+allruns[irun]+$
            '/spec1d_unfluxed/scra.*fits',count=nfile)
          outfiles = atlaspath+allruns[irun]+'/spec1d/'+$
            file_basename(repstr(allfiles,'scra.','fscra.'))
          psfile = atlaspath+allruns[irun]+'/spec1d/qa_'+$
            allruns[irun]+'.ps'
          im_plotconfig, 8, pos, psfile=psfile
          for ii = 0, nfile-1 do begin
             long_fluxcal, allfiles[ii], sensfunc=sensfuncfile, $
               outfil=outfiles[ii]
; make the QAplot
             flux = mrdfits(outfiles[ii],0,hdr)
             wave = mrdfits(outfiles[ii],2)
             galaxy = sxpar(hdr,'OBJECT')
             frame = repstr(file_basename(outfiles[ii]),'.fits','')
             frame = strmid(frame,strpos(frame,'.')+1)
             djs_plot, wave, flux, xsty=3, ysty=3, ps=10, $
               position=pos, xtitle='Wavelength (\AA)', $
               ytitle='Flux (10^{-17} '+flam_units()+')'
             legend, [galaxy,frame], /right, /top, box=0
          endfor
          im_plotconfig, psfile=psfile, /psclose, /gzip
       endif
       
; --------------------------------------------------
; copy and rename the final 1D spectra to a central clearinghouse
       if keyword_set(final_spec1d) then begin
          dwave = 2.75D         ; output dispersion
          trim = 10             ; remove ten pixels on either end
          singsfile = atlaspath+allruns[irun]+'/singslist.txt'
          if file_test(singsfile,/reg) then begin
             readcol, singsfile, infile, outfile, galaxy, phot, $
               /silent, format='A,A,A,A', comment='#'
             nobj = n_elements(infile)
; grab the redshifts for these objects from NED
             allgalaxy = repstr(galaxy,'_','')
             ugalaxy = allgalaxy[uniq(allgalaxy,sort(allgalaxy))]
             ned_webget_basic, ugalaxy, ned
; now loop on each spectrum             
             for ii = 0, nobj-1 do begin
                infile1 = atlaspath+allruns[irun]+'/spec1d/'+infile[ii]
                outfile1 = spec1dpath+strlowcase(strtrim(outfile[ii],2))+'.fits'
                if file_test(infile1,/reg) then begin
                   flux = mrdfits(infile1,0,hdr)
                   ferr = mrdfits(infile1,1)
                   wave = mrdfits(infile1,2)
                   npix = n_elements(wave)
                   flux = 1E-17*flux[trim:npix-trim-1]
                   ferr = 1E-17*ferr[trim:npix-trim-1]
                   wave = wave[trim:npix-trim-1]
; rebin linearly in wavelength
                   newwave = dindgen((max(wave)-min(wave))/dwave+1)*$
                     dwave+min(wave)
                   newflux = rebin_spectrum(flux,wave,newwave)
                   newvar = rebin_spectrum(ferr^2,wave,newwave)
                   newferr = sqrt(newvar*(newvar gt 0.0)) ; enforce positivity
; build the final header; also grab the redshift from NED
                   sxaddpar, hdr, 'CRVAL1', min(wave), ' wavelength at CRPIX1'
                   sxaddpar, hdr, 'CRPIX1', 1.0D, ' reference pixel number'
                   sxaddpar, hdr, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'
                   sxaddpar, hdr, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
                   sxaddpar, hdr, 'CTYPE1', 'LINEAR', ' projection type'
                   sxaddpar, hdr, 'PHOTFLAG', phot[ii], ' photometric flag [Y/N]'
                   sxaddpar, hdr, 'GALAXY', galaxy[ii], ' galaxy name'
                   thisned = where(repstr(galaxy[ii],'_','') eq ugalaxy)
                   sxaddpar, hdr, 'ZNED', float(ned[thisned].z), ' NED redshift'
                   plot, newwave, newflux, ps=10, xsty=3, ysty=3
; write out                   
                   if file_test(outfile1,/reg) then begin
                      splog, 'Output file '+file_basename(outfile1)+' exists!'
;                     outfile1 = 
                   endif
                   mwrfits, newflux, outfile1, hdr, /create
                   mwrfits, newferr, outfile1, hdr
                   spawn, 'gzip -f '+outfile1, /sh
                endif else begin
                   splog, 'File '+infile[ii]+' not yet reduced!'
                endelse
             endfor 
             stop
          endif
       endif 

       stop          
       
    endfor                      ; run

    
    

return
end
