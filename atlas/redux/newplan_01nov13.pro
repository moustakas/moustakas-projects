function newplan_01nov13, planfile, hdr=hdr, sings=sings, $
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
      "xtrim '0,0'",$
      'prof_nsigma 50',$
;     'nolocal 1',$
      'stddir Std']

    scidir = where(strmatch(hdr,'*scidir*'))
    if keyword_set(sings) then hdr[scidir] = 'scidir SINGS' else $
      hdr[scidir] = 'scidir ATLAS'

; the twiflats are mis-identified as standards
    these = where(strmatch(old.target,'*skyflat*'))
    new[these].flavor = 'twiflat'

; remove the "nuclear" spectrum of NGC1560, because there is no
; well-defined nucleus
    rej = where($
      (strtrim(new.filename,2) eq 'a.3746.fits.gz') or $
      (strtrim(new.filename,2) eq 'a.3747.fits.gz'),comp=keep)
    new = new[keep]

; remove the 2.5" and non-CALSPEC standard stars
    new = long_find_calspec(new,datapath)
    keep = where((strtrim(new.flavor,2) ne 'std') or $
      ((strtrim(new.flavor,2) eq 'std') and $
      strmatch(new.maskname,'*4.5*') and $
      (strmatch(new.starname,'*...*') eq 0)))
    new = new[keep]
    
; get the aperture info and return
    apinfo = atlas_get_apertures(new,sings=sings,atlas=atlas)

return, new
end
