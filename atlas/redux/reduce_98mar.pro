pro reduce_98mar, plan=plan, reduce=reduce, sensfunc=sensfunc, $
  coadd=coadd, fluxcal=fluxcal, clobber=clobber

    planfile = 'plan_98mar.par'

; --------------------------------------------------
    if keyword_set(plan) then begin
       long_plan, '*.fits.gz', 'Raw/', planfile=planfile
       old = yanny_readone(planfile,hdr=hdr)
       new = old

; turn off optimal extraction and allow for extended sources 
       hdr = hdr[where(strcompress(hdr,/remove) ne '')]
       maxobj = where(strmatch(hdr,'*maxobj*'))
       hdr[maxobj] = 'maxobj 5'
       hdr = [hdr,$
         'slitthresh 100.0',$
         'prof_nsigma 20',$
         'nolocal 1']

; the twiflats are mis-identified as standards
       these = where(strmatch(old.filename,'a.010[5-7].fits.gz'))
       new[these].flavor = 'twiflat'

; fix MASKNAME for the afternoon calibration data
       these = where(strmatch(old.filename,'a.010[1-7].fits.gz'))
       new[these].maskname = '2.5'

; do not reduce the 2.5" standards; also drop non-calspec standards 
       rej = where((strmatch(new.flavor,'*std*') and $
         strmatch(new.maskname,'*2.5*')) or $
         strmatch(new.filename,'a.0110.fits*'),comp=keep)
       new = new[keep]

; write out the new plan file
       yanny_write, planfile, ptr_new(new), hdr=hdr, /align
       return
    endif else plan = yanny_readone(planfile)

; --------------------------------------------------
; reduce    
    if keyword_set(reduce) then begin
       long_reduce, planfile, /verbose, clobber=clobber
    endif

; --------------------------------------------------
; make the sensitivity function
    sensfuncfile = 'sensfunc_98mar.fits'
    if keyword_set(sensfunc) then begin
       
       std = where((strtrim(plan.flavor,2) eq 'std'),nstd)
       stdfiles = 'Science/std-'+strtrim(plan[std].filename,2)
       info = iforage(stdfiles)
       pad = where(strmid(info.ra,0,1) eq ' ',npad)
       if (npad ne 0) then info[pad].ra = '0'+strtrim(info[pad].ra,2)

       allstd = rsex(getenv('XIDL_DIR')+'/Spec/Longslit/'+$
         'calib/standards/calspec/calspec_info.txt')
       spherematch, 15.0*im_hms2dec(allstd.ra), im_hms2dec(allstd.dec), $
         15.0*im_hms2dec(info.ra), im_hms2dec(info.dec), 1.0, $
         m1, m2, max=0
       std_names = repstr(allstd[m1].file,'.fits.gz')

       sens = im_long_sensfunc(stdfiles[m2],sensfuncfile,$
         std_name=std_names,/msk_balm)
    endif

; --------------------------------------------------
; need: crsplits_ files
    date = '98mar21'
    if keyword_set(coadd) then begin
       coaddfile = 'crsplits_'+date+'.txt'
       coaddtxt = djs_readlines(coaddfile)
       coaddtxt = coaddtxt[where((strcompress(coaddtxt,/remove) ne '') and $
         strmatch(coaddtxt,'*#*') eq 0)]
       for ii = 0, n_elements(coaddtxt)-1 do begin
          infiles = strtrim(strsplit(coaddtxt[ii],' ',/extract),2)
          infiles = repstr(infiles,'ra.','sci-a.')+'.gz' ; convert ispec file name
          nobs = n_elements(infiles)
          outfile = strmid(infiles[0],0,strpos(infiles[0],'.fits.gz'))+$
            '_'+string(nobs,format='(I0)')+'.fits'
          outfile = repstr(outfile,'sci-','scr') ; hack!
          long_coadd, 'Science/'+infiles, 1, outfil='Fspec/'+outfile
       endfor
    endif
       
; --------------------------------------------------
; flux-calibrate
    if keyword_set(fluxcal) then begin
       allfiles = file_search('Fspec/scra.*fits',count=nfile)
       outfiles = repstr(allfiles,'scra.','fscra.')
       for ii = 0, nfile-1 do begin
          long_fluxcal, allfiles[ii], sensfunc=sensfuncfile, $
            outfil=outfiles[ii]
       endfor
    endif

stop    

return
end

;    long_slitmask, 'Raw/a.0105.fits.gz', 'test.fits', $
;      biasfile='superbias-a.0101.fits', tset_slits=tset, $
;      xstart=5, xend=120
