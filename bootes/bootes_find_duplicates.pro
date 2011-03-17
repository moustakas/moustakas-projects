pro bootes_find_duplicates
; jm10jan27ucsd - simple code to test M. Brown's FLAG_DUPLICATE switch
; jm10oct12ucsd - rewritten to find and write out unique observations

    common bootes, cat1, ingroup, mult, first, next
    
    year = '2010b' ; '2009b'
    iopath = getenv('DATA_DIR')+'/data/bootes/'+year+'/'
    ibandfile = iopath+'bootes_I.fits.gz'
    outfile = repstr(ibandfile,'.fits.gz','.uniq.fits')

    if (n_elements(cat1) eq 0) then begin
       splog, 'Reading '+ibandfile
       cat1 = mrdfits(ibandfile,1)
    endif
    ngal = n_elements(cat1)
    
; use spheregroup to remove duplicate observations
    cat = struct_addtags(cat1,replicate({ndup: 0, isbest: 0},ngal))

    if (n_elements(ingroup) eq 0L) then begin
       ingroup = spheregroup(cat.alpha_j2000,cat.delta_j2000,$
         1D/3600D,multgroup=mult,firstgroup=first,nextgroup=next)
    endif
    ugroup = ingroup[uniq(ingroup,sort(ingroup))] ; unique groups
    nugroup = n_elements(ugroup)
    splog, 'Parsing '+string(nugroup,format='(I0)')+' unique groups'
;   finalindx = lindgen(ngal)

    t0 = systime(1)
    for igroup = 0L, nugroup-1L do begin
       if ((igroup mod 1000) eq 0) then splog, 'group '+$
         string(igroup,format='(I6.6)')

       match = where((ingroup eq ugroup[igroup]),nmatch)
       if (nmatch eq 1L) then cat[match].isbest = 1 ; just one
       if (nmatch gt 1) then begin
         struct_print, /no_head, struct_trimtags(cat[match],$
           select=['alpha_j2000','delta_j2000','mag_auto',$
           'magerr_auto','flag_duplicate',$
           'flag_subfield','image','class_star']), ddigit=11
         print

; check the various flags; from M. Brown:          

; SEGFLAGS indicates the number of pixels assigned to other objects
; within the aperture. Unfortunately it is in units of pixels rather
; than square arcsec. A big value indicates that lots of pixels within
; the aperture have been assigned to other objects and while the code
; attempts to correct for this, it does an imperfect job.

; IMAFLAGS provides the minimum value of the exposure or weight map
; within the aperture. For the SDWFS and NDWFS data it is the minimum
; number of individual exposures contributing to a pixel. Anything >=3
; should be good although you could apply more conservative
; criteria. For the NEWFIRM data it is the weight map multiplied by
; 1000, and anything >=100 seems to be OK.

; FLAG_SUBFIELD indicates that an object is within the nominal
; subfield boundaries. Objects with FLAG_SUBFIELD=0 may be near the
; edges of images and be of poorer quality. I tend to use this flag in
; my object selection.

; FLAG_DUPLICATE=1 indicates that the object is almost certainly a
; duplicate detection of an object already in the
; catalogue. Duplicates arise because of the overlap of NDWFS
; subfields, and I use this flag for my object selection.

; FLAG_DEEP=1 indicates that deep optical data is available at the
; relevant location. This gets around the problem that bright star
; halos and the broken up components of galaxies are not normally
; flagged by IRAF (because they are not saturated etc). The catch with
; FLAG_DEEP is all bright objects (e.g., I<17) will tend to have
; FLAG_DEEP=0 so you may need to ignore the FLAG_DEEP keyword for the
; brightest objects. I tend to use flag deep when dealing with faint
; object samples (and Mark Brodwin has found it very useful too.)

; first check the subfield flag
         subf = where(cat[match].subfield eq 1,nsubf)
         if (nsubf eq 1) then begin
            cat[match[subf]].isbest = 1
            continue
         endif else begin

stop            
         endelse
         
         snr = max(cat[match].flux_psf/cat[match].fluxerr_psf,maxsnr)
         cat[match].ndup = nmatch
         cat[match[maxsnr]].isbest = 1

         if (igroup gt 50000) then stop
       endif
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60.0

stop    
    
; write out    
    im_mwrfits, cat, outfile, /clobber

return
end
