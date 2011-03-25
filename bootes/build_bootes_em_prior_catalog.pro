;+
; NAME:
;   build_bootes_em_prior_catalog
; PURPOSE:
;   Build the optical prior photometric catalog for use with the
;   GALEX/EM code. 
; CALLING SEQUENCE:
;
; INPUTS: 
;
; KEYWORDS:
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Jul 26, UCSD
;-

pro build_bootes_em_prior_catalog, combine=combine, gr=gr, clobber=clobber

    if (n_elements(gr) eq 0) then gr = 'gr6'
    empath = getenv('IM_ARCHIVE_DIR')+'/bootes/galex/'
    photodir = getenv('IM_ARCHIVE_DIR')+'/bootes/2010b/'

; --------------------------------------------------
; combine the Bw- and I-band catalogs, remove duplicates, and write
; out a combined catalog that forms the basis for constructing the
; prior catalog
    if keyword_set(combine) then begin
       iband1 = mrdfits(photodir+'bootes_I.fits.gz',1)
       bband1 = mrdfits(photodir+'bootes_Bw.fits.gz',1)

; remove duplicates
       keep = where((iband1.flag_subfield eq 1) and (iband1.flag_duplicate eq 0),ngal)
       iband = iband1[keep]
       bband = bband1[keep]

; use the psf photometry
       out = im_struct_trimtags(iband,select=['alpha_j2000','delta_j2000',$
         'mag_psf','class_star'],newtags=['ra','dec','imag','class_star'])
       out = struct_addtags(temporary(out),im_struct_trimtags(bband,$
         select='mag_psf',newtags='bmag'))

       crap = where((out.imag le 0.0) or (out.imag gt 90.0),ncrap)
       out[crap].imag = iband[crap].mag_auto

       im_mwrfits, out, empath+'priors/bootes_combined.fits', /clobber
    endif

; --------------------------------------------------
; now build the prior catalog    
    survey = 'ngpdws'
    priorcat1 = {ident: 0L, ra: 0.0D, dec: 0.0D, $
      u: 0.0D, r: 0.0D, star: 0, magflag: 0}

; read the GR6 tiles that overlap the PRIMUS survey regions and just
; keep the DIS ones for now
    tileinfo = read_galex_tileinfo()
    tileinfo = tileinfo[where(strtrim(tileinfo.survey,2) eq 'DIS')]
    struct_print, tileinfo
    
    ratag = 'ra'
    dectag = 'dec'
    umagtag = 'bmag'
    rmagtag = 'imag'
    photofile = 'bootes_combined.fits.gz'

; output file names       
    priorfile = empath+'priors/'+survey+'.prior.fits'
    psfile = repstr(priorfile,'.fits','.ps')
    infofile = repstr(priorfile,'.prior.fits','.galextiles.dat')
    if file_test(priorfile+'.gz') and (keyword_set(clobber) eq 0) then begin
       splog, 'Output file '+priorfile+'.gz exists; use /CLOBBER'
       continue
    endif
       
; read the optical catalog
    splog, 'Reading '+photodir+photofile
    cat = mrdfits(photodir+photofile,1)

    cat = struct_addtags(cat,replicate({ident: 0L},n_elements(cat)))
    cat.ident = lindgen(n_elements(cat))
    raindx = tag_indx(cat,ratag)
    decindx = tag_indx(cat,dectag)
    rmagindx = tag_indx(cat,rmagtag)
    umagindx = tag_indx(cat,umagtag)

; spherematch against the galex tiles and only keep the area of
; interest
    spherematch, cat.(raindx), cat.(decindx), tileinfo.tile_ra, $
      tileinfo.tile_dec, 0.65, m1, m2, maxmatch=0
    m2 = m2[uniq(m2,sort(m2))]
    splog, 'Found '+strtrim(n_elements(m2),2)+' GALEX tiles'
    keep = where($
      (cat.(raindx) gt min(tileinfo[m2].tile_ra)-0.65) and $
      (cat.(raindx) lt max(tileinfo[m2].tile_ra)+0.65) and $
      (cat.(decindx) lt max(tileinfo[m2].tile_dec)+0.65) and $
      (cat.(decindx) gt min(tileinfo[m2].tile_dec)-0.65),nobj)
    
    cat = cat[keep]
    
; apply some basic quality cuts and find stars       
    good = lindgen(n_elements(cat))
    isstar = (cat[good].imag lt 19.0) and (cat[good].class_star gt 0.95)
    ngood = n_elements(good)
    outcat = cat[good]

    priorcat = replicate(priorcat1,ngood)
    priorcat.ident = outcat.ident
    priorcat.ra = outcat.(raindx)
    priorcat.dec = outcat.(decindx)       
    priorcat.u = outcat.(umagindx)
    priorcat.r = outcat.(rmagindx)
    priorcat.star = isstar

; following discussions with S. Arnouts, if an object does *not* have
; a "u-band" detection but *does* have an "r-band" detection, then
; replace the u-band magnitude with the r-band magnitude and add a
; flag indicating that this has been done
    rep = where(((priorcat.u le 0.0) or (priorcat.u gt 90.0)) and $
      ((priorcat.r gt 0.0) and (priorcat.r lt 90.0)),nrep)
    if (nrep ne 0L) then begin
       splog, 'Replacing '+strtrim(nrep,2)+' magnitudes'
       priorcat[rep].magflag = 1
       priorcat[rep].u = priorcat[rep].r
    endif

; write out and then make the QAplot       
    splog, 'Writing '+priorfile
    mwrfits, priorcat, priorfile, /create
    struct_print, tileinfo[m2], file=infofile

    bin = ngood/20
    im_plotconfig, 0, pos, psfile=psfile, width=6.5, $
      height=6.5, xmargin=[1.5,0.5], ymargin=[0.6,1.1]
    djs_plot, priorcat.ra, priorcat.dec, position=pos, $
      psym=3, symsize=0.1, xsty=3, ysty=3, bin=bin, $
      xtitle='\alpha_{J2000} (degree)', $
      ytitle='\delta_{J2000} (degree)', color='grey'
    tt = tileinfo[m2]
;      tt = tileinfo[sort(tileinfo.tile_ra)]
    for itile = 0, n_elements(tt)-1 do begin
       if odd(itile) then offset = -0.1 else offset = 0.1
       tvcircle, 0.65, tt[itile].tile_ra, tt[itile].tile_dec, $
         color=djs_icolor('blue'), /data, thick=5.0
       xyouts, tt[itile].tile_ra, tt[itile].tile_dec+offset, $
         strtrim(tt[itile].tile_name,2), align=0.5, charsize=1.6
    endfor
    im_plotconfig, /psclose, /gzip, psfile=psfile
    
return
end
