pro unpack_gdds, doplot=doplot, wfits=wfits
; jm05feb08uofa
; unpack the GDDS data
    
    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    msun = +5.47                ; the Sun's absolute B-band magnitude
    
    rootpath = gdds_path(/raw)
    outpath = gdds_path(/spec1d)
    zpath = gdds_path(/analysis)
    
; read the redshift data

    zfile = 'GDDS-z.cat.sex'
    splog, 'Reading '+zpath+zfile+'.'
    zinfo = rsex(zpath+zfile)

; only include relevant objects

    keep = where((zinfo.zconf ge 2L) and (zinfo.zconf le 9L) and (zinfo.agn eq 0L) and $
      (zinfo.coll ge 0L) and (zinfo.coll le 2L) and (zinfo.z_obj ge 0.4) and (zinfo.z_obj le 1.0),nobj)
    zinfo = zinfo[keep]
;   struct_print, struct_trimtags(zinfo[keep],select=['ID','Z_OBJ','ZCONF','COLL','AGN','OII','OIII','CLASS'])

    galaxy = strtrim(zinfo.id,2)
    objlist = galaxy+'.txt'
    path = strmid(galaxy,0,4)
    
    pagemaker, nx=1, ny=2, ymargin=[0.2,1.0], yspace=0.0, $
      xmargin=[1.3,0.2], position=position, /normal

; loop on each object

    for iobj = 4L, nobj-1L do begin
;   for iobj = 0L, nobj-1L do begin

       if (not keyword_set(wfits)) then $
         print, format='("Processing object ",I0,"/",I0,".",A1,$)', iobj+1, nobj, string(13b)

       datapath = rootpath+path[iobj]+'/'
       zmatch = zinfo[iobj]

; does this object have an error spectrum?  if not then generate an
; error spectrum according to the mean S/N
       
         line = djs_readilines(datapath+objlist[iobj],indx=0L)
         if strmatch(line,'*NOTE: THIS IS AN EXPORTED *EXTERNAL* SPECTRUM*') then begin

            readcol, datapath+objlist[iobj], objwave, objflux, $
              format='D,D', comment='#', /silent
            good = where(objflux ne 0.0,ngood)
            objwave = objwave[good]
            objflux = objflux[good]

; compute the median S/N according to the statistics in the continuum
; at the central wavelength of the spectrum            
            
            snrwave = djs_mean(objwave)
            get_element, objwave, snrwave+30*[-1,+1], indx
            specsnr = median(objflux[indx[0]:indx[1]])/djsig(objflux[indx[0]:indx[1]],sigrej=3.0) ; median S/N
            
; next generate the error spectrum according to the square root of the
; counts in the spectrum            

            get_element, objwave, snrwave, normindx
            objferr = objflux/sqrt((objflux/objflux[normindx])^2)/specsnr
            
         endif else begin

; read the data and error spectrum

            readcol, datapath+objlist[iobj], objwave, objflux, objferr, $
              format='D,D,D', comment='#', /silent

         endelse

       npix = n_elements(objwave)
       
; regrid the data onto a constant wavelength array

       minwave = 3727*(zmatch.z_obj+1)-100.0 > min(objwave)
       maxwave = 5007*(zmatch.z_obj+1)+100.0 < max(objwave)
       dwave = (maxwave - minwave) / (npix - 1.0)

       newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave

       good = where(objferr lt 1E10)
       objflux = interpol(objflux[good],objwave[good],newwave)
       objferr = sqrt(interpol(objferr[good]^2,objwave[good],newwave))
       objwave = newwave

; interpolate over negative pixels

       neg = where(objflux le 0.0,nneg,comp=pos,ncomp=npos)

; correct for foreground Galactic extinction       

       ra = im_hms2dec(zmatch.ra)*15.0
       dec = im_hms2dec(zmatch.dec)
       
       glactc, ra, dec, 2000.0, gl, gb, 1, /degree
       ebv_mw = dust_getval(gl,gb,/interp)

       kl = k_lambda(objwave,/odonnell,R_V=3.1)

       objflux = objflux * 10^(0.4*kl*ebv_mw)
       objferr = objferr * 10^(0.4*kl*ebv_mw)

       if keyword_set(doplot) then begin

          xrange = minmax(objwave)

          ploterror, objwave, objflux, objferr, ps=10, xsty=3, ysty=3, position=position[*,0], /normal, $
;         djs_plot, objwave, objflux, ps=10, xsty=3, ysty=3, position=pos[*,0], /normal, $
            thick=2.0, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
            ytitle='Object [Counts]', charsize=1.5, charthick=2.0, xrange=xrange
          legend, galaxy[iobj], /left, /top, box=0, charsize=1.5, charthick=2.0
          if (nneg ne 0L) then djs_oplot, objwave[neg], objflux[neg], ps=4, color='red'
          
          djs_plot, objwave, objflux/objferr, ps=10, xsty=3, ysty=3, position=position[*,1], /normal, $
            /noerase, xthick=2.0, ythick=2, ytitle='S/N', charsize=1.5, charthick=2.0, $
            xtitle='Observed Wavelength [\AA]', xrange=xrange

          cc = get_kbrd(1)
          
       endif

       mkhdr, header, objflux, /extend
       sxdelpar, header, 'COMMENT'

       sxaddpar, header, 'RA', zmatch.ra
       sxaddpar, header, 'DEC', zmatch.dec
       sxaddpar, header, 'GALAXY', galaxy[iobj]
       sxaddpar, header, 'Z', float(zmatch.z_obj), ' redshift'
       sxaddpar, header, 'Z_ERR', float(0.0), ' redshift error'
       sxaddpar, header, 'CRVAL1', minwave, ' wavelength at CRPIX1'
       sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
       sxaddpar, header, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
       sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
       sxaddpar, header, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'

       if keyword_set(wfits) and (zmatch.z_obj gt 0.0) and (zmatch.z_obj lt 1.0) then begin

          outfile = galaxy[iobj]+'.fits'

          splog, 'Writing '+outpath+outfile+'.'
          mwrfits, objflux, outpath+outfile, header, /create
          mwrfits, objferr, outpath+outfile
          spawn, ['gzip -f '+outpath+outfile], /sh
          
       endif

    endfor 

stop
    
return
end
