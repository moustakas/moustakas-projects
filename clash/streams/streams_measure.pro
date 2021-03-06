;+
; NAME:
;   dmeasure_atlas
; PURPOSE:
;   Run measurements on atlas results
; CALLING SEQUENCE:
;   dmeasure_atlas, subdir
; COMMENTS:
;   Assumes dimage-style detection outputs
; REVISION HISTORY:
;   31-July-2010
;-

pro streams_measure, base, psffiles=psffiles, bulgedisk=bulgedisk, $
  mgefit=mgefit, noclobber=noclobber

    sub='atlases'
    postfix=''

; default to use base name same as directory name
    if n_elements(base) eq 0 then begin
       spawn, 'pwd', cwd
       base=(file_basename(cwd))[0]
    endif

; read in pset
    pset= gz_mrdfits(base+'-pset.fits',1)
    spawn, 'pwd', subdir, /sh
    subdir=subdir[0]
    
    ref= (pset.ref)[0]
    nbands= n_elements(pset.imfiles)

; hack to find pixel scale
    pixscales= fltarr(nbands)
    for i=0L, nbands-1L do begin
       hdr= gz_headfits(subdir+'/'+strtrim(pset.imfiles[i],2))
       nx= long(sxpar(hdr, 'NAXIS1'))
       ny= long(sxpar(hdr, 'NAXIS2'))
       ntest=10L
       xyad, hdr, nx/2L, ny/2L, ra1, dec1
       xyad, hdr, nx/2L+ntest, ny/2L, ra2, dec2
       cirrange,ra1
       cirrange,ra2
       spherematch, ra1, dec1, ra2,dec2, 360., m1, m2, d12
       pixscales[i]=d12/float(ntest)
    endfor
    scales= pixscales[ref]/pixscales

; psf file names
    if n_elements(psffiles) eq 0 then begin
       psffiles= strarr(nbands)
       for i=0L, nbands-1L do begin
          imfile= strtrim(string(pset.imfiles[i]),2)
          psffiles[i]=(strmid(imfile, 0, strlen(imfile)-8)+'-bpsf.fits')
       endfor
    endif

; loop on each parent
    pcat=gz_mrdfits(base+'-pcat.fits',1)

    splog, 'Hack!!'
;   for iparent = 266, 275 do begin
;   for iparent = 13, 20 do begin
;   for iparent = 54L, n_elements(pcat)-1 do begin
    for iparent = 0L, n_elements(pcat)-1 do begin
       splog, 'Parent ', iparent
    
       pstr = strtrim(string(iparent),2)
       
       mfile=subdir+'/'+sub+'/'+pstr+'/'+base+'-'+pstr+ $
         '-measure'+postfix+'.fits'
       sfile=subdir+'/'+sub+'/'+pstr+'/'+base+'-'+pstr+ $
         '-sersic'+postfix+'.fits'
       
       if(keyword_set(noclobber) eq 0 OR $
         gz_file_test(mfile) eq 0) then begin
          
          acat=gz_mrdfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-acat-'+ $
            pstr+'.fits',1)
          nkeep=0
          
          if(n_tags(acat) gt 0) then ikeep=where(acat.good gt 0 and acat.type eq 0,nkeep)

; loop on each object in the catalog
          if nkeep gt 0 then begin
             acat=acat[ikeep]

             for ik = 0, nkeep-1 do begin
                aid=acat[ik].aid
                astr=strtrim(string(aid),2)
             
                rimage=gz_mrdfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-'+ pstr+'-atlas-'+astr+'.fits', ref, rhdr)
                rinvvar=gz_mrdfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-'+ $
                  'ivar-'+pstr+'.fits', ref)
                rinvvar=rinvvar>0.

                if ((keyword_set(noclobber) eq 0 OR gz_file_test(mfile) eq 0) AND $
                  keyword_set(nomeasure) eq 0) then begin
                
                   adxy, rhdr, acat[ik].racen, acat[ik].deccen, xcen, ycen
                   psf= gz_mrdfits(psffiles[ref])

                   dmeasure, rimage, rinvvar, xcen=xcen, ycen=ycen, $
                     measure=r_measure
                   r_sersic=0

; keep the "sky" as a free parameter to allow for the ICL!                    
                   if keyword_set(mgefit) then begin
                      stop                      
                   endif else begin
                      if keyword_set(bulgedisk) then begin
                         dsersic2, rimage, rinvvar, xcen=r_measure.xcen, ycen=r_measure.ycen, $
                           sersic=r_sersic, /fixcen, model=refmodel, psf=psf, $
                           bulge=bulge, disk=disk, fixsky=1
                      endif else begin
;if ik eq 1 then stop
                         dsersic, rimage, rinvvar, xcen=r_measure.xcen, ycen=r_measure.ycen, $
                           sersic=r_sersic, /fixcen, model=refmodel, psf=psf, fixsky=0
                      endelse
                   endelse
; this is a bug in Blanton's code, I think
;                  dsersic, rimage, rinvvar, xcen=xcen, ycen=ycen, sersic=r_sersic, $
;                    /fixcen, /fixsky, model=refmodel, psf=psf

                   if ik eq 0 then outmodel = fltarr([size(rimage,/dim),nbands])

; comment this out                   
;                  outhdr= rhdr
;                  sxdelpar, outhdr, 'XTENSION'
;                  sxdelpar, outhdr, 'PCOUNT'
;                  sxdelpar, outhdr, 'GCOUNT'
;                  mwrfits, float(refmodel), sfile, outhdr, /create
                   
                   xyad, rhdr, r_measure.xcen, r_measure.ycen, racen, deccen
                   cirrange, racen
                
                   help,/st,r_measure
;                  print, xcen, ycen, r_measure.xcen, r_measure.ycen, $
;                    r_sersic.xcen, r_sersic.ycen
                   
                   mall= {racen:racen, $
                     deccen:deccen, $
                     xcen:r_measure.xcen, $
                     ycen:r_measure.ycen, $
                     nprof:bytarr(nbands), $
                     profmean:fltarr(nbands, 15), $
                     profmean_ivar:fltarr(nbands, 15), $
                     profradius:r_measure.profradius, $
                     qstokes:fltarr(nbands, 15), $
                     ustokes:fltarr(nbands, 15), $
                     bastokes:fltarr(nbands, 15), $
                     phistokes:fltarr(nbands, 15), $
                     petroflux:fltarr(nbands), $
                     petroflux_ivar:fltarr(nbands), $
;                    fiberflux:fltarr(nbands), $
;                    fiberflux_ivar:fltarr(nbands), $
                     petrorad:r_measure.petrorad, $
                     petror50:r_measure.petror50, $
                     petror90:r_measure.petror90, $
                     ba50:r_measure.ba50, $
                     phi50:r_measure.phi50, $
                     ba90:r_measure.ba90, $
                     phi90:r_measure.phi90, $
                     chi2: fltarr(nbands), $
                     sky: fltarr(nbands), $
                     sersicflux:fltarr(nbands), $
                     sersicflux_ivar:fltarr(nbands), $
                     sersic_r50:r_sersic.sersicr50, $
                     sersic_n:r_sersic.sersicn, $
                     sersic_ba:r_sersic.axisratio, $
                     sersic_phi:r_sersic.orientation, $
                     asymmetry:fltarr(nbands), $
                     clumpy:fltarr(nbands), $
                     dflags:lonarr(nbands), $
                     aid:aid}

                   for iband=0L, nbands-1L do begin
                      image=gz_mrdfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-'+ $
                        pstr+'-atlas-'+astr+'.fits', iband, hdr,/silent)
                      invvar=gz_mrdfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-'+ $
                        'ivar-'+pstr+'.fits', iband,/silent)
                      invvar=invvar>0.
                      
                      psf = gz_mrdfits(psffiles[iband], /silent)
                   
                      adxy, hdr, racen, deccen, xcen, ycen
                      dmeasure, image, invvar, xcen=xcen, ycen=ycen, /fixcen, $
                        measure=tmp_measure, cpetrorad= mall.petrorad*scales[iband], $
                        faper=7.57576*scales[iband]
                      curr_sersic= r_sersic
                      curr_sersic.xcen= xcen
                      curr_sersic.ycen= ycen
                      curr_sersic.sersicr50= r_sersic.sersicr50*scales[iband]

                      delvarx, thismodel
                      if keyword_set(mgefit) then begin
                         
                      endif else begin
                         if keyword_set(bulgedisk) then begin
                            dsersic2, image, invvar, xcen=xcen, ycen=ycen, $
                              sersic=curr_sersic, /fixcen, /onlyflux, model=thismodel, $
                              psf=psf, bulge=bulge, disk=disk, fixsky=1
                         endif else begin
                            dsersic, image, invvar, xcen=xcen, ycen=ycen, $
                              sersic=curr_sersic, /onlyflux, /fixcen, fixsky=0, $
                              psf=psf, model=thismodel
                         endelse
                      endelse
                      thismodel = thismodel - curr_sersic.sky ; subtract the sky!
                      outmodel[*,*,iband] = outmodel[*,*,iband] + thismodel
                      
                      mall.nprof[iband]= tmp_measure.nprof
                      mall.profmean[iband,*]= tmp_measure.profmean
                      mall.profmean_ivar[iband,*]= $
                        tmp_measure.profmean_ivar
                      mall.qstokes[iband,*]= tmp_measure.qstokes
                      mall.ustokes[iband,*]= tmp_measure.ustokes
                      mall.bastokes[iband,*]= tmp_measure.bastokes
                      mall.phistokes[iband,*]= tmp_measure.phistokes
                      mall.petroflux[iband]= tmp_measure.petroflux
                      mall.petroflux_ivar[iband]= $
                        tmp_measure.petroflux_ivar
;                     mall.fiberflux[iband]= tmp_measure.fiberflux
;                     mall.fiberflux_ivar[iband]= tmp_measure.fiberflux_ivar
                      mall.chi2[iband]= curr_sersic.chisquared/(curr_sersic.ndof+$
                        (curr_sersic.ndof eq 0))*(curr_sersic.ndof gt 0)
                      mall.sky[iband]= curr_sersic.sky
                      mall.sersicflux[iband]= curr_sersic.sersicflux
                      mall.sersicflux_ivar[iband]= curr_sersic.sersicflux_ivar
                      mall.asymmetry[iband]= tmp_measure.asymmetry
                      mall.clumpy[iband]= tmp_measure.clumpy
                      mall.dflags[iband]= tmp_measure.dflags
                   endfor ; close loop on each band
;                  dhdr= dimage_hdr()
;                  mwrfits, mall, mfile, dhdr, /create
;                  spawn, 'gzip -f '+mfile, /sh
                endif 
                if ik eq 0 then outmall = mall else outmall = [outmall,mall]
             endfor             ; close loop on each object
; write out fitting results
             dhdr= dimage_hdr()
             mwrfits, outmall, mfile, dhdr, /create
             spawn, 'gzip -f '+mfile, /sh
; write out Sersic model fits
             for iband=0L, nbands-1L do begin
                outhdr=gz_headfits(subdir+'/'+sub+'/'+pstr+'/'+base+'-'+ $
                  pstr+'-atlas-'+astr+'.fits',ext=iband)
                sfile=subdir+'/'+sub+'/'+pstr+'/'+base+'-'+pstr+ $
                  '-sersic'+postfix+'.fits'
                mwrfits, reform(float(outmodel[*,*,iband])), $
                  sfile, outhdr, create=iband eq 0
             endfor
             spawn, 'gzip -f '+sfile, /sh
          endif 
       endif 
    endfor ; close parent loop

return
end
