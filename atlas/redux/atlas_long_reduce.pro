;+
; NAME:
;   atlas_long_reduce
;
; PURPOSE:
;
;   Main program for the Low-redux pipeline.  This set of algorithms
;   runs mainly as a black box.
;
; CALLING SEQUENCE:
;  long_reduce, planfile, /clobber, /NOZAP, /NOFLEX, /NOHELIO 
;
; INPUTS:
;  planfile  -- File created by long_plan which guides the reduction
;               process
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;               pixels)
;  PROF_NSIGMA= -- Extend the region to fit a profile by hand 
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
; LINELIST=  -- File for arc calibration lines
; REID_FILE=  -- File for cross-correlation [Replaces the default].
;                Note you should set BIN_RATIO when doing this.
; BIN_RATIO=  -- Ratio of binning in the data over the binning of the archived arc
; ARC_INTER=  -- Fiddle with the wavelengths interactively:  
;                1 -- Run x_identify with the archive template 
;                2 -- Run x_identify without it
;
; OUTPUTS:
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;
; BUGS:
;   
; REVISION HISTORY:
;   11-Mar-2005  Written by JH + SB
;   J. Moustakas, 2009-Dec-19, UCSD - various edits
;-  
;-----------------------------------------------------------------------------
PRO atlas_long_reduce, planfile, clobber = clobber, verbose = verbose $
                 , NOFLEX = noflex, NOZAP = NOZAP, NOHELIO = NOHELIO $
                 , ONLYSCI = ONLYSCI, HAND_X = HAND_X, HAND_Y = HAND_Y $
                 , HAND_FWHM = HAND_FWHM, LINELIST = linelist $
                 , REID_FILE= REID_FILE $
                 , USAGE = USAGE, FILESTD = FILESTD1, CHK = CHK $
                 , TRCCHK = TRCCHK, NOSHIFT = NOSHIFT, STD = STD $
                 , ISLIT = ISLIT, PROFILE = PROFILE, _EXTRA = extra $
                 , SKYTRACE=SKYTRACE, PROF_NSIGMA=prof_nsigma $
                 , BIN_RATIO=bin_ratio, ARC_INTER=arc_inter, $
  justcalib=justcalib, juststd=juststd, justsci=justsci, $
  calibclobber=calibclobber, novac=novac, maxobj=maxobj, ccdproc=ccdproc, $ ; jm09dec19ucsd
  badpixfile=badpixfile

    
  if  KEYWORD_SET(USAGE) THEN BEGIN
      print,'Syntax - ' + $
        'long_reduce, planfile, /CLOBBER, /VERBOSE, /NOFLEX, /NOHELIO [v1.0]' 
      return
  endif 

if (NOT keyword_set(planfile)) then planfile = findfile('plan*.par')

   ;----------
   ; If multiple plan files exist, then call this script recursively
   ; for each such plan file.

   if planfile[0] EQ '' then begin
      print, 'could not find plan file'
      print, 'try running long_plan'
      return
   endif

   if (n_elements(planfile) GT 1) then begin
      for i=0L, n_elements(planfile)-1L do begin
          splog, '------------------------'
          splog, 'Reducing ' + planfile[i]
          splog, '------------------------'
          long_reduce, planfile[i], clobber = clobber, verbose = verbose $
                       , NOFLEX = noflex, NOZAP = NOZAP, NOHELIO = NOHELIO $
                       , ONLYSCI = ONLYSCI, HAND_X = HAND_X, HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM, LINELIST = linelist $
                       , USAGE = USAGE, FILESTD = FILESTD1, CHK = CHK $
                       , NOSHIFT = NOSHIFT, STD = STD $
                       , ISLIT = ISLIT, PROFILE = PROFILE, _EXTRA = extra
      endfor
      return
  endif

   ;----------
   ; Read the plan file

   planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
   if (NOT keyword_set(planstr)) then begin
       splog, 'Empty plan file ', planfile
       return
   endif

   ;; Truncate the gz stuff
   nfil = n_elements(planstr)
   for qq=0L,nfil-1 do begin
       slen = strlen(planstr[qq].filename)
       if strmid(planstr[qq].filename,slen-3) EQ '.gz' then $
         planstr[qq].filename = strmid(planstr[qq].filename,0,slen-3)
   endfor
   logfile = yanny_par(planhdr, 'logfile')
   plotfile = yanny_par(planhdr, 'plotfile')
   indir = yanny_par(planhdr, 'indir')
   tempdir = yanny_par(planhdr, 'tempdir')
   scidir  = yanny_par(planhdr, 'scidir')
   stddir  = yanny_par(planhdr, 'stddir')
   preprocdir = yanny_par(planhdr, 'preprocdir')
;   write_flats = yanny_par(planhdr, 'write_flats')
;   subsample = yanny_par(planhdr, 'subample')
   minslit1 = yanny_par(planhdr, 'minslit')
   slitthresh1 =  yanny_par(planhdr, 'slitthresh')
   reduxthresh = yanny_par(planhdr, 'reduxthresh')
   sig_thresh = yanny_par(planhdr, 'sig_thresh')
   nolocal = yanny_par(planhdr, 'nolocal')
   prof_nsigma = yanny_par(planhdr, 'prof_nsigma')
   slity1_1 = yanny_par(planhdr, 'slity1')
   slity2_1 = yanny_par(planhdr, 'slity2')
   maxflex1 = yanny_par(planhdr, 'maxflex')
   maxgood = yanny_par(planhdr, 'maxgood')
;   ksize_1  = yanny_par(planhdr, 'slitksize')

   box_rad1 = yanny_par(planhdr, 'box_rad')
   if (n_elements(maxobj) eq 0) then $ ; jm09dec19ucsd
     maxobj = yanny_par(planhdr, 'maxobj')
   xtrim = yanny_par(planhdr, 'xtrim') ; jm09dec19ucsd

   IF NOT KEYWORD_SET(STD) THEN STD = yanny_par(planhdr, 'std') 
   IF KEYWORD_SET(slity1_1) THEN slity1 = long(slity1_1)
   IF KEYWORD_SET(slity2_1) THEN slity2 = long(slity2_1)
;   IF KEYWORD_SET(ksize_1) THEN ksize = long(ksize_1)
   IF KEYWORD_SET(box_rad1) THEN box_rad = long(box_rad1)
   IF KEYWORD_SET(minslit1) THEN minslit = long(minslit1)
   IF KEYWORD_SET(maxflex1) THEN maxflex = long(maxflex1)

;   IF KEYWORD_SET(slitthresh1) THEN slitthresh = long(slitthresh1)

   ;----------
   ; Create science dir
   IF keyword_set(scidir) THEN spawn, '\mkdir -p '+scidir
   IF keyword_set(stddir) THEN spawn, '\mkdir -p '+stddir
;  IF keyword_set(preprocdir) THEN spawn, '\mkdir -p '+preprocdir ; jm09dec28ucsd

   ;----------
   ; Open log file
;  if keyword_set(justcalib) then $ ; jm09dec19ucsd
;    logfile = repstr(logfile,'.log','.calib.log') 
   if keyword_set(juststd) then $ ; jm09dec19ucsd
     logfile = repstr(logfile,'.log','.std.log') 
   if (keyword_set(logfile)) then begin
       splog, filename = logfile
       splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'idlutils version ' + idlutils_version()
   splog, 'Longslit version ' + longslit_version()
   
   plotfile = 0
;  if keyword_set(onlycalib) then $ ; jm09dec19ucsd
;    plotfile = repstr(plotfile,'.ps','.calib.ps') 
   if keyword_set(juststd) then $ ; jm09dec19ucsd
     plotfile = repstr(plotfile,'.ps','.std.ps') 
   if (keyword_set(plotfile)) then begin
       thisfile = findfile(plotfile, count = ct)
       IF (ct EQ 0 OR KEYWORD_SET(CLOBBER)) THEN BEGIN
           splog, 'Plot file ' + plotfile
           dfpsplot, plotfile, /color
       ENDIF ELSE BEGIN
           cpbackup, plotfile
           splog, 'Plot file already exists. Creating backup'
           splog, 'Plot file ' + plotfile
           dfpsplot, plotfile, /color
       ENDELSE
   ENDIF

   ;----------
   ; Loop over each INSTRUMENT (e.g., CCD)

   ccd_list = planstr.instrument
   ccd_list = ccd_list[uniq(ccd_list, sort(ccd_list))]
   nccd = n_elements(ccd_list)

   for iccd=0L, nccd-1L do begin
      indx = where(planstr.instrument EQ ccd_list[iccd])

      ;;-----------------------------
      ;; Make a superbias if possible
      ;;-----------------------------

      ii = where(planstr[indx].flavor EQ 'bias', nbias)
      if (nbias GT 0) then begin
          ibias = indx[ii]
          superbiasfile = 'superbias-' + planstr[ibias[0]].filename
          thisfile = findfile(superbiasfile, count = ct)
          if (ct EQ 0 OR keyword_set(calibclobber)) then begin
              splog, 'Generating superbias for INSTRUMENT=', ccd_list[iccd]
              long_superbias $
                , djs_filepath(planstr[ibias].filename $
                               , root_dir = indir) $
                , superbiasfile, verbose = verbose
          endif else begin
              splog, 'Do not overwrite existing superbias ', thisfile
          endelse
      endif else begin
          superbiasfile = ''
          splog, 'No input biases for superbias for INSTRUMENT=', $
                 ccd_list[iccd]
      endelse
  
      ; Loop over each GRATING+MASK+WAVE for this INSTRUMENT

      mask_list = planstr[indx].grating 
; + planstr[indx].wave + planstr[indx].maskname $
      mask_list = mask_list[uniq(mask_list, sort(mask_list))]
      nmask = n_elements(mask_list)
      
      for imask = 0L, nmask-1L do begin
          jndx = indx[ where(planstr[indx].grating  EQ mask_list[imask]) ]
          
          qboth = planstr[jndx].flavor EQ 'bothflat'
          qtwi = (planstr[jndx].flavor EQ 'twiflat') OR qboth
          qpix = (planstr[jndx].flavor EQ 'domeflat' OR $
                   planstr[jndx].flavor EQ 'iflat') OR qboth
;          iflat = where(qtwi OR qdome)
          itwi  = where(qtwi, ntwi)
          ipix  = where(qpix, npix)
          iboth = WHERE(qtwi OR qpix, nboth)
         ;;---------------------------
         ;; Find the slits on the mask and create slitmask file 
         ;;---------------------------
          if (npix GT 0 OR ntwi GT 0) then begin
              ;; Use twiflats if possible, otherwise domeflats.
              ;; (Use only the first such file).
              if (ntwi GT 0) then ithis = jndx[itwi[0]] $
              else ithis = jndx[ipix[0]] 
              slitfile = 'slits-' + planstr[ithis].filename
              thisfile = findfile(slitfile, count = ct)
              if (ct EQ 0 OR keyword_set(calibclobber)) then begin
                  splog, 'Generating slits for INSTRUMENT=', ccd_list[iccd], $
                         ' GRATING+MASK+WAVE=', mask_list[imask]
                  IF KEYWORD_SET(SLITTHRESH1) THEN BEGIN
                      slitthresh = double(slitthresh1) 
                      nfind = 0
                  ENDIF ELSE IF $
                    strmatch(planstr[ithis].MASKNAME, '*long*') THEN BEGIN
                      slitthresh = 0.15D
                      ;; JXP -- Allow for LRISb  5/30/08
                      if strmatch(strtrim(ccd_list[iccd],2),'LRISBLUE') then $
                        nfind = 2 else begin
                          ;; Red side (allow for new chips)
                          hdr = xheadfits(djs_filepath(planstr[ithis].filename,$
                                                      root_dir=indir))
                          if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then $
                            nfind = 1 else nfind = 2
                      endelse
                  ENDIF ELSE BEGIN
                      slitthresh = 0.02
                      nfind = 0
                  ENDELSE
                  long_slitmask, djs_filepath(planstr[ithis].filename, $
                                              root_dir = indir) $
                                 , slitfile, minslit = minslit $
                                 , peakthresh = slitthresh $
                                 , y1 = slity1, y2 = slity2 $
                                 , ksize = ksize, nfind = nfind $
                                 , biasfile = superbiasfile, verbose = verbose, $
                    xtrim=xtrim ; jm09dec19ucsd
              endif else begin
                  splog, 'Do not overwrite existing slitmask file ', $
                         thisfile
              endelse
          endif else begin
              slitfile = ''
              splog, 'No input flats for slitmask for INSTRUMENT=', $
                     ccd_list[iccd], ' GRATING+MASK+WAVE=', mask_list[imask]
          endelse
     
         ;---------------------------
         ; Make a wavelength solution
         ;---------------------------

; jm09dec19ucsd - reduce all the arc lamps
         ii = where(planstr[jndx].flavor EQ 'arc', narc)

         if (narc GT 0) then begin
             iarc = jndx[ii]
             for jarc = 0, narc-1 do begin
                wavefile = 'wave-' + planstr[iarc[jarc]].filename ; jm09dec18ucsd
                if strpos(wavefile,'.fits') LT 0 then wavefile=wavefile+'.fits'
                thisfile = findfile(wavefile, count = ct)
                if (ct EQ 0 OR keyword_set(calibclobber)) then begin
                   splog, 'Generating wavelengths for INSTRUMENT=' $
                     , ccd_list[iccd], $
                     ' GRATING+MASK+WAVE=', mask_list[imask]
                   long_wavesolve $
                     , djs_filepath(planstr[iarc[jarc]].filename $ ; jm09dec18ucsd
                     , root_dir = indir), wavefile $
                     , slitfile = slitfile, biasfile = superbiasfile $
                     , verbose = verbose, LINELIST = linelist, CHK = CHK $
                     , REID_FILE=reid_file, BIN_RATIO=bin_ratio, $
                     ARC_INTER=arc_inter
                endif else begin
                   splog, 'Do not overwrite existing wavelength file ', thisfile
                endelse
             endfor
          endif else begin
             wavefile = ''
             splog, 'No input arcs for wavelengths for INSTRUMENT=', $
               ccd_list[iccd], ' GRATING+MASK+WAVE=', mask_list[imask]
          endelse

         ;;------------------------------------
         ;; Make  pixel and illumination flats 
         ;;------------------------------------
         IF (npix GT 0) OR (ntwi GT 0) and $
           (not strmatch(planstr[0].instrument,'GMOS-NX')) AND $
           (not strmatch(planstr[0].instrument,'GMOS-SX')) THEN BEGIN 
             IF npix GT 0 THEN BEGIN
                 pixflatfile = 'pixflat-' + planstr[jndx[ipix[0]]].filename
                 thispixflatfile = findfile(pixflatfile, count = pixct)
             ENDIF ELSE pixct = 0
             IF ntwi GT 0 THEN BEGIN
                 illumflatfile = 'illumflat-' + planstr[jndx[itwi[0]]].filename
                 thisillumflatfile = findfile(illumflatfile, count = illumct)
             ENDIF ELSE aillumct = 0
             if (pixct EQ 0 OR illumct EQ 0 OR keyword_set(calibclobber)) then begin
                 ;; Higher order for MMT illumination effect
                 ;; F strmatch(ccd_list[0], 'MMT*') OR $
                 ;;  strmatch(ccd_list[0], 'DEIMOS*') THEN npoly = 7L
                 splog, 'Generating pixel flat for INSTRUMENT=' $
                        , ccd_list[iccd], ' GRATING+MASK+WAVE=' $
                        , mask_list[imask]
                 if keyword_set(calibclobber) then begin ; jm09nov06ucsd
                    pixct = 0
                    illumct = 0
                 endif
                 IF illumct EQ 0 AND pixct EQ 0 THEN BEGIN
                     infiles =  djs_filepath(planstr[jndx[iboth]].filename $
                                             , root_dir = indir)
;                    use_illum = qtwi
;                    use_pixel = qpix
                     use_illum = qtwi[iboth] 
                     use_pixel = qpix[iboth] 
                 ENDIF ELSE IF illumct EQ 1 AND pixct EQ 0 THEN BEGIN
                     infiles =  djs_filepath(planstr[jndx[ipix]].filename $
                                             , root_dir = indir)
                     use_illum = lonarr(npix)
                     use_pixel = lonarr(npix) + 1L
                     splog, 'Do not overwrite existing illumination flat ' $
                            , thisillumflatfile
                 ENDIF ELSE IF illumct EQ 0 AND pixct EQ 1 THEN BEGIN
                     if ntwi NE 0 then begin 
                         infiles =  djs_filepath(planstr[jndx[itwi]].filename $
                                                 , root_dir = indir) 
                         use_illum = lonarr(ntwi) + 1L
                         use_pixel = lonarr(ntwi)
                     endif else begin
                        infiles =  $
                           djs_filepath(planstr[jndx[ipix[0]]].filename $
                                        , root_dir = indir) 
                         use_illum = 0L
                         use_pixel = 0L
                      endelse
                     splog, 'Do not overwrite existing pixel flats ' $
                            , thispixflatfile
                 ENDIF

                 long_superflat, infiles, pixflatfile, illumflatfile $
                                 , slitfile = slitfile $
                                 , wavefile = wavefile $
                                 , biasfile = superbiasfile $
                                 , verbose = verbose $
                                 , tempdir = tempdir $
                                 , use_illum = use_illum $
                                 , _EXTRA = extra $
                                 , use_pixel = use_pixel $
                                 , npoly = npoly, CHK = CHK
             ENDIF ELSE BEGIN
                 splog, 'Do not overwrite existing pixel flats ' $
                        , thispixflatfile
                 splog, 'Do not overwrite existing illumination flat ' $
                        , thisillumflatfile
             ENDELSE
         endif else begin
             pixflatfile = ''
             splog, 'No input pixel flats or illum flats for INSTRUMENT=', $
                    ccd_list[iccd], ' GRATING+MASK=', mask_list[imask]
         endelse

;          IF strmatch(ccd_list[0], 'LRISBLUE') THEN BEGIN
;             IF mask_list[imask] eq '300/5000' THEN BEGIN
               
;                splog, 'Correcting LRISBLUE 300/5000 pixflat with pixels of unity below 4000 Angstroms.'
;                print, 'Correcting LRISBLUE 300/5000 pixflat with pixels of unity below 4000 Angstroms.'

;                ;Read in the wavelength solution file already created.
;                wavefile_fix = 'wave-' + planstr[iarc[0]].filename
;                if strpos(wavefile_fix,'.fits') LT 0 then wavefile_fix=wavefile_fix+'.fits'
;                wavefile_fix = findfile(wavefile_fix, count = ct_fix)
;                IF ct_fix GT 0 then wavefile_fix_im = mrdfits(wavefile_fix[0], 0) else stop

;                ;Read in the pixel flat file already created.
;                pixflatfile_fix = 'pixflat-' + planstr[jndx[ipix[0]]].filename
;                pixflatfile_fix = findfile(pixflatfile_fix, count = pixct_fix)
;                if pixct_fix gt 0 then pixflatfile_fix_im = mrdfits(pixflatfile_fix[0], 0, headfix) else stop
               
;                scattered = where(wavefile_fix_im lt 4000. and wavefile_fix_im ne 0., ct_scatt)
;                if ct_scatt gt 0 then pixflatfile_fix_im[scattered]=1. else stop
               
;                print, 'Deleting old pixflat file for LRISBLUE 300/5000 setup.'
;                spawn, 'rm ' + pixflatfile[0]
;                mwrfits, pixflatfile_fix_im, pixflatfile_fix[0], headfix

;             ENDIF
;          ENDIF

; don't reduce the science data         
         if keyword_set(justcalib) then continue ; jm09dec19ucsd

         planstr1 = planstr
; just reduce the standards
         if keyword_set(juststd) then begin ; jm09dec19ucsd
            planstr1 = planstr
            keep = where(strtrim(planstr1.flavor,2) eq 'std',nkeep)
            if (nkeep eq 0) then begin
               splog, 'No standards observed'
               return
            endif
            planstr = planstr1[keep]
            indx = where(planstr.instrument EQ ccd_list[iccd])
            jndx = indx[ where(planstr[indx].grating  EQ mask_list[imask]) ]
            nozap = 1 ; very important!!!!
         endif

; just reduce the science data
         if keyword_set(justsci) then begin ; jm09dec19ucsd
; choose a standard star as a trace crutch
            istd = where(planstr[jndx].flavor EQ 'std', nstd)
            IF nstd GT 0 THEN BEGIN
               stdfile = 'std-' + planstr[jndx[istd[0]]].filename
               thisstdfile = findfile(stddir+'/' + stdfile + '*', count = ct1)
;              thisstdfile = findfile(scidir+'/' + stdfile + '*', count = ct1)
               IF ct1 GT 0 THEN filestd = thisstdfile[0]
            endif

; now crop the plan to just include the science objects
            planstr1 = planstr
            keep = where(strtrim(planstr1.flavor,2) eq 'science',nkeep)
            if (nkeep eq 0) then begin
               splog, 'No science objects observed'
               return
            endif
            planstr = planstr1[keep]
            indx = where(planstr.instrument EQ ccd_list[iccd])
            jndx = indx[ where(planstr[indx].grating  EQ mask_list[imask]) ]
         endif

         ;-----------------------------------
         ; Finally, reduce each science image
         ;-----------------------------------

         ii = where(planstr[jndx].flavor EQ 'science' OR $
                    planstr[jndx].flavor EQ 'std', nsci)
         for isci = 1, 1 do begin
;        for isci = 0L, nsci-1 do begin
             j = jndx[ii[isci]]
             IF KEYWORD_SET(ONLYSCI) THEN $
               IF strmatch(ONLYSCI, planstr[j].FILENAME) NE 1 THEN CONTINUE
             if planstr[j].flavor EQ 'science' then begin
                if keyword_set(ccdproc) then $ ; jm09dec20ucsd
                  scifile = djs_filepath('ccd-'+planstr[j].filename,root_dir = preprocdir) else $
                    scifile = djs_filepath('sci-' + planstr[j].filename,root_dir = scidir)
             endif else begin
                scifile = djs_filepath('std-' + planstr[j].filename, root_dir = stddir)
;               scifile = djs_filepath('std-' + planstr[j].filename, root_dir = scidir)
             endelse
             ipos = strpos(scifile, '.fits')
             if ipos LT 0 then scifile = scifile+'.fits'
;             IF KEYWORD_SET(PROFILE) THEN $
;               profile_filename = 'profile-' + planstr[j].filename $
;             ELSE PROFILE_FILENAME= 0
             thisfile = findfile(scifile+'*', count = ct)
             if (ct EQ 0 OR keyword_set(clobber)) THEN BEGIN
                 splog, 'Reducing science frame ', prelog = planstr[j].filename
                 IF KEYWORD_SET(FILESTD1) THEN filestd = filestd1 $
                 ELSE IF planstr[j].flavor EQ 'science' then begin
                     istd = where(planstr[jndx].flavor EQ 'std', nstd)
                     IF nstd GT 0 THEN BEGIN
                         stdfile = 'std-' + planstr[jndx[istd[0]]].filename
                         thisstdfile = findfile(stddir+'/' + stdfile + '*', count = ct1)
;                        thisstdfile = findfile(scidir+'/' + stdfile + '*', count = ct1)
                         IF ct1 GT 0 THEN filestd = thisstdfile[0] 
                     ENDIF 
                 ENDIF 
; choose the wavelength calibration map nearest in time to the
; observations; jm09dec28ucsd
                 allarc = where(planstr1.flavor EQ 'arc',nallarc)
                 windx = atlas_choose_wavefile(indir+planstr[j].filename,$
                   indir+planstr1[allarc].filename)
                 wavefile = 'wave-'+planstr1[allarc[windx]].filename
                 splog, 'Choosing WAVEFILE '+wavefile

; basic processing; jm09dec28ucsd
                 if keyword_set(ccdproc) then begin
                    atlas_long_ccdproc, djs_filepath(planstr[j].filename $
                      , root_dir = indir), scifile $
                      , slitfile = slitfile, wavefile = wavefile $
                      , biasfile = superbiasfile $
                      , pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile $
                      , verbose = verbose $
                      , maxobj = maxobj, box_rad = box_rad $
                      , reduxthresh = reduxthresh $
                      , sig_thresh = sig_thresh $
                      , noflex = noflex, NOZAP = NOZAP $
                      , maxflex = maxflex , MAXGOOD=maxgood $
                      , NOHELIO = NOHELIO $
;                                   , profile_filename = profile_filename $
                      , HAND_X = HAND_X, HAND_Y = HAND_Y $
                      , HAND_FWHM = HAND_FWHM $
                      , PROF_NSIGMA=PROF_NSIGMA $
                      , FILESTD = FILESTD $
                      , STD = (planstr[j].flavor EQ 'std') $
                      , ISLIT = ISLIT, CHK = CHK $
                      , TRCCHK = TRCCHK, NOLOCAL=nolocal $
                      , NOSHIFT = NOSHIFT, SKYTRACE = SKYTRACE $
                      , _EXTRA = extra, $
                      novac=novac, badpixfile=badpixfile ; jm09dec19ucsd
                    continue
                 endif

;; deal with the hand apertures                 
;                 if tag_exist(planstr,'hand_x') then begin
;                    notzero = where(planstr[j].hand_x gt 0.0,nnotzero)
;                    if (nnotzero ne 0) then begin
;                       hand_x = planstr[j].hand_x[notzero]
;                       hand_y = planstr[j].hand_y[notzero]
;                       hand_fwhm = planstr[j].hand_fwhm[notzero]
;                       box_rad = planstr[j].box_rad[notzero]
;                    endif
;                 endif

;                hand_x = 60
;                hand_y = 600
;                hand_fwhm = 50.0
                 
                 atlas_long_reduce_work, djs_filepath(planstr[j].filename $
                   , root_dir = indir), scifile $
                                   , slitfile = slitfile, wavefile = wavefile $
                                   , biasfile = superbiasfile $
                                   , pixflatfile = pixflatfile $
                                   , illumflatfile = illumflatfile $
                                   , verbose = verbose $
                                   , maxobj = maxobj, box_rad = box_rad $
                                   , reduxthresh = reduxthresh $
                                   , sig_thresh = sig_thresh $
                                   , noflex = noflex, NOZAP = NOZAP $
                                   , maxflex = maxflex , MAXGOOD=maxgood $
                                   , NOHELIO = NOHELIO $
;                                   , profile_filename = profile_filename $
                                   , HAND_X = HAND_X, HAND_Y = HAND_Y $
                                   , HAND_FWHM = HAND_FWHM $
                                   , PROF_NSIGMA=PROF_NSIGMA $
                                   , FILESTD = FILESTD $
                                   , STD = (planstr[j].flavor EQ 'std') $
                                   , ISLIT = ISLIT, CHK = CHK $
                                   , TRCCHK = TRCCHK, NOLOCAL=nolocal $
                                   , NOSHIFT = NOSHIFT, SKYTRACE = SKYTRACE $
                                   , _EXTRA = extra, $
                   novac=novac, badpixfile=badpixfile ; jm09dec19ucsd
                 splog, prelog = ''
             endif else begin
                 splog, 'Do not overwrite existing science frame ', scifile
             endelse
          endfor 
   endfor                       ; End loop over GRATING+MASK
   endfor ; End loop over INSTRUMENT

   if (keyword_set(plotfile) and !d.name eq 'PS') then begin
       dfpsclose
   endif
   
   splog, /close
   x_psclose
   
   return
end
;------------------------------------------------------------------------------
