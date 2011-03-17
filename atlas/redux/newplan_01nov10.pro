function newplan_01nov10, planfile, hdr=hdr, sings=sings, $
  atlas=atlas, apinfo=apinfo
; jm09dec17ucsd - edit the PLAN file

    datapath = file_dirname(planfile)+'/'
    old = yanny_readone(planfile,hdr=hdr)
    new = old
    
; turn off optimal extraction and allow for extended sources 
    hdr = hdr[where(strcompress(hdr,/remove) ne '')]

    maxobj = where(strmatch(hdr,'*maxobj*'))
    hdr[maxobj] = 'maxobj 5'

    reduxthresh = where(strmatch(hdr,'*reduxthresh*'))
    hdr[reduxthresh] = 'reduxthresh 0.1'
    
    hdr = [hdr,$
;     'slitthresh 100.0',$
;     "xtrim '5,10'",$
      "xtrim '0,0'",$
      'prof_nsigma 100',$
      'nolocal 1',$
      'stddir Std']
;     'preprocdir Preproc']

    scidir = where(strmatch(hdr,'*scidir*'))
    if keyword_set(sings) then hdr[scidir] = 'scidir SINGS' else $
      hdr[scidir] = 'scidir ATLAS'

; the twiflats are mis-identified as standards
    these = where(strmatch(old.target,'*skyflat*'))
    new[these].flavor = 'twiflat'

; scan error on a.3426 (NGC 584 drift 55): remove from analysis; also
; reject a.3428 (NGC 584 nuc), which was not at the parallactic angle
    rej = where($
      (strtrim(new.filename,2) eq 'a.3426.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3428.fits.gz'),comp=keep)
    new = new[keep]

; remove the 2.5" and non-CALSPEC standard stars
    new = long_find_calspec(new,datapath)
    keep = where((strtrim(new.flavor,2) ne 'std') or $
      ((strtrim(new.flavor,2) eq 'std') and $
      strmatch(new.maskname,'*4.5*') and $
      (strmatch(new.starname,'*...*') eq 0)))
    new = new[keep]

; ---------------------------------------------------------------------------
; SINGS/ATLAS galaxies that need to be stitched are handled
; separately, so write out a dedicated plan file here
    stitch = where($
      (strtrim(new.filename,2) eq 'a.3418.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3419.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3420.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3421.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3438.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3439.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3440.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3441.fits.gz'),comp=keep)

    if keyword_set(sings) then begin
       stitch_plan = new[cmset_op(stitch,'or',where(strtrim(new.flavor,2) ne 'science'))]
       stitch_planfile = repstr(planfile,'.par','.stitch.par')
       splog, 'Writing '+stitch_planfile
       yanny_write, stitch_planfile, ptr_new(stitch_plan), hdr=hdr, /align
    endif

    new = new[keep]
; ---------------------------------------------------------------------------

; get the aperture info and return
    apinfo = atlas_get_apertures(new,sings=sings,atlas=atlas)

return, new
end
