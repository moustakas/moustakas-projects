pro streams_bcg
; jm13sep04siena - do the preliminary BCG subtraction

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    struct_print, sample
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]
    rmaxkpc = 60D     ; [kpc]

; wrap on each cluster    
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       skypath = streams_path(/skysub)+cluster+'/'
       outpath = streams_path(/nobcg)+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       rmax = ceil(rmaxkpc/arcsec2kpc/pixscale)         ; [pixels]

; read the skyinfo structure to get the filters
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster[ic]+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       reffilt = where(short eq 'f160w') ; reference filter
       nfilt = n_elements(skyinfo)
       niceprint, skyinfo.band, skyinfo.sigma, skyinfo.sblimit
    
; read and fit the galaxy in the reference image          
       scales = fltarr(nfilt)+1 ; this could be generalized

       psffiles = file_search(streams_path(/psfs)+short+'-bpsf.fits*')
       psf = gz_mrdfits(psffiles[reffilt])

       imfile = skypath+cluster+'-'+short[reffilt]+'.fits.gz'
       splog, 'Reading '+imfile
       rimage = mrdfits(imfile,0,rhdr,/silent)
       rinvvar = mrdfits(imfile,1,/silent)

       bcgra1 = 15D*hms2dec(sample[ic].ra)
       bcgdec1 = hms2dec(sample[ic].dec)
       adxy, rhdr, bcgra1, bcgdec1, xcen, ycen
       
       r_sersic = 0
       dmeasure, rimage, rinvvar, xcen=xcen, ycen=ycen, measure=r_measure
       dsersic, rimage, rinvvar, xcen=r_measure.xcen, ycen=r_measure.ycen, $
         sersic=r_sersic, /fixcen, model=refmodel, psf=psf, /fixsky
;      dsersic2, rimage, rinvvar, xcen=r_measure.xcen, ycen=r_measure.ycen, $
;        sersic=r_sersic, /fixcen, model=refmodel, psf=psf, bulge=bulge, $
;        disk=disk, /fixsky

       xyad, rhdr, r_measure.xcen, r_measure.ycen, racen, deccen
       cirrange, racen
          
; now loop through each band and scale the Sersic fit             
       for ib = nfilt-1, 0, -1 do begin
;      for ib = 0, nfilt-1 do begin
          imfile = skypath+cluster+'-'+short[ib]+'.fits.gz'
          splog, 'Reading '+imfile
          image = gz_mrdfits(imfile,0,hdr,/silent)
          invvar = gz_mrdfits(imfile,1,ivarhdr,/silent)
          invvar = invvar>0.0

          psf = gz_mrdfits(psffiles[ib],/silent)

          adxy, hdr, racen, deccen, xcen, ycen
          curr_sersic = r_sersic
          curr_sersic.xcen = xcen
          curr_sersic.ycen = ycen
          curr_sersic.sersicr50 = r_sersic.sersicr50*scales[ib]
             
          dsersic, image, invvar, xcen=xcen, ycen=ycen, sersic=curr_sersic, $
            /onlyflux, /fixcen, /fixsky, psf=psf, model=model
;         dsersic2, image, invvar, xcen=xcen, ycen=ycen, sersic=curr_sersic, $
;           /nofit, /fixcen, /fixsky, psf=psf, model=model, bulge=bulge, disk=disk
             
;; get a RMAXKPC by RMAXKPC cutout centered on the BCG
;          extast, hdr, astr
;          ad2xy, clash[this].ra_bcg, clash[this].dec_bcg, astr, xcen, ycen
;          xcen = fix(xcen)
;          ycen = fix(ycen)
;          hextract, image, hdr, cutimage, cuthdr, xcen-rmax, $
;            xcen+rmax, ycen-rmax, ycen+rmax, /silent
;          hextract, invvar, ivarhdr, cutinvvar, cutivarhdr, xcen-rmax, $
;            xcen+rmax, ycen-rmax, ycen+rmax, /silent
;          badpixels = where(cutinvvar le 0,nbad)
;          if nbad eq 0L then delvarx, badpixels
;
;          mge1 = streams_mge(cutimage,badpixels=badpixels,$
;            pixscale=pixscale,model=cutmodel,twist=twist)
;          help, mge1, /str

; rebuild the original image with the BCG subtracted and write out 
;         model = image*0.0
;         embed_stamp, model, cutmodel, xcen-rmax, ycen-rmax
;         atv, image-model, /bl

          outfile = outpath+cluster[ic]+'-'+short[ib]+'.fits'
          splog, 'Writing '+file_basename(outfile)

          mwrfits, image-model, outfile, hdr, /create
          mwrfits, invvar, outfile, ivarhdr, /silent
          spawn, 'gzip -f '+outfile
       endfor                   ; close filter loop
    endfor                      ; close cluster loop

return
end
    
