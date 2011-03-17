pro unpack_gdds, doplot=doplot, wfits=wfits
; jm05feb08uofa
; unpack the GDDS data
    
    red, h100=0.7, omega_0=0.3, omega_lambda=0.7
    msun = +5.47                ; the Sun's absolute B-band magnitude
    
    rootpath = gdds_path(/raw)
    outpath = gdds_path(/spec1d)
    zpath = gdds_path(/analysis)
    
;   subdirs = ['SA02','SA12','SA22']
;   subdirs = ['SA02','SA12','SA15','SA22']
;   nsubdirs = n_elements(subdirs)
    
; read the redshift data

    zfile = 'GDDS-z.cat.sex'
    splog, 'Reading '+zpath+zfile+'.'
    zinfo = rsex(zpath+zfile)

; only include relevant objects

    keep = where((zinfo.zconf ge 2L) and (zinfo.zconf le 9L) and (zinfo.agn eq 0L) and $
      (zinfo.coll ge 0L) and (zinfo.coll le 2L) and (zinfo.z_obj ge 0.4) and (zinfo.z_obj le 1.0))
    zinfo = zinfo[keep]
;   struct_print, struct_trimtags(zinfo[keep],select=['ID','Z_OBJ','ZCONF','COLL','AGN','OII','OIII','CLASS'])

    galaxy = strtrim(zinfo.id,2)
    objlist = galaxy+'.txt'
    datapath = strmid(galaxy,0,4)
    
; loop on each subdirectory    

    pagemaker, nx=1, ny=2, ymargin=[0.2,1.0], yspace=0.0, $
      xmargin=[1.3,0.2], position=position, /normal

;   for idir = 3L, nsubdirs-1L do begin
    for idir = 0L, nsubdirs-1L do begin
    
; construct the object list

       datapath = rootpath+subdirs[idir]+'/'
       splog, 'Pushing into '+datapath
    
       pushd, datapath
       objlist = file_search('*.txt',count=nobj)

       popd

; exclude objects with no error spectrum

       keep = lindgen(nobj)

       delvarx, nosigma
       for iobj = 0L, nobj-1L do begin

          line = djs_readilines(datapath+objlist[iobj],indx=0L)
          if strmatch(line,'*NOTE: THIS IS AN EXPORTED *EXTERNAL* SPECTRUM*') then begin
             if (n_elements(nosigma) eq 0L) then nosigma = iobj else nosigma = [nosigma,iobj]
             
          endif

       endfor

       if (n_elements(nosigma) ne 0L) then begin

          remove, nosigma, keep
          objlist = objlist[keep]
          nobj = n_elements(objlist)

       endif

; loop on each object

;      for iobj = 40, nobj-1L do begin
       for iobj = 0L, nobj-1L do begin

          if (not keyword_set(wfits)) then $
            print, format='("Processing object ",I0,"/",I0,".",A1,$)', iobj+1, nobj, string(13b)

; extract the root object name

          galaxy = strtrim(strupcase(repstr(objlist[iobj],'.txt','')),2)

; match this object to the redshift data file

          match = where(strmatch(zinfo.id,galaxy,/fold) eq 1B,nmatch)

          if (nmatch eq 0L) then begin
             print
             splog, '**** No REDSHIFT match for '+objlist[iobj]+' ****'
             print
          endif
          if (nmatch gt 1L) then begin
             print
             splog, '**** Multiple matches for '+objlist[iobj]+' ****'
             print
          endif

          if (nmatch eq 1L) then begin

             zmatch = zinfo[match]

; read the data and error spectrum

             readcol, datapath+objlist[iobj], objwave, objflux, objferr, $
               format='D,D,D', comment='#', /silent
             npix = n_elements(objwave)

; regrid the data onto a constant wavelength array

             case galaxy of
                'SA02-1741': minwave = 5400.0
                'SA02-1778': minwave = 5700.0
;               'SA12-5224': minwave = 5200.0
                'SA12-6301': minwave = 5200.0
                else: minwave = min(objwave)
             endcase

             case galaxy of
                'SA02-1878': maxwave = 8700.0
                'SA02-1933': maxwave = 7400.0
                'SA12-6232': maxwave = 8400.0
                'SA22-2476': maxwave = 8500.0
                else: maxwave = max(objwave) ; 9900.0 < max(objwave)
             endcase
             
             dwave = (maxwave - minwave) / (npix - 1.0)

             newwave = findgen((maxwave-minwave)/dwave+1)*dwave+minwave

             good = where(objferr lt 1E10)
             objflux = interpol(objflux[good],objwave[good],newwave)
             objferr = sqrt(interpol(objferr[good]^2,objwave[good],newwave))
             objwave = newwave

; add a constant offset equal to the most negative pixel value, plus
; 10%, but preserve the S/N

;            oldflux = objflux
;            oldferr = objferr

;            objflux = objflux + abs(min(objflux))*1.1
;            objferr = objferr * abs(oldflux/objflux)
             
;            snr_median_old = strtrim(string(median(oldflux/oldferr),format='(F12.2)'),2)
;            snr_median = strtrim(string(median(objflux/objferr),format='(F12.2)'),2)
;            print, snr_median_old+' '+snr_median
             
; interpolate over negative pixels

             neg = where(objflux le 0.0,nneg,comp=pos,ncomp=npos)
;            if (npos eq 0L) then message, 'Problem here!'
;            
;            if (nneg ne 0L) then begin
;
;               while (nneg ne 0L) do begin
;
;                  objflux = interpol(objflux[pos],objwave[pos],objwave)
;                  objferr = sqrt(interpol(objferr[pos]^2,objwave[pos],objwave))
;                  neg = where(objflux le 0.0,nneg)
;
;               endwhile
;               
;            endif
             
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
;               djs_plot, objwave, objflux, ps=10, xsty=3, ysty=3, position=pos[*,0], /normal, $
                  thick=2.0, xthick=2.0, ythick=2.0, xtickname=replicate(' ',10), $
                  ytitle='Object [Counts]', charsize=1.5, charthick=2.0, xrange=xrange
                legend, galaxy, /left, /top, box=0, charsize=1.5, charthick=2.0
                if (nneg ne 0L) then djs_oplot, objwave[neg], objflux[neg], ps=4, color='red'
                
                djs_plot, objwave, objflux/objferr, ps=10, xsty=3, ysty=3, position=position[*,1], /normal, $
                  /noerase, xthick=2.0, ythick=2, ytitle='S/N', charsize=1.5, charthick=2.0, $
                  xtitle='Observed Wavelength ['+angstrom()+']', xrange=xrange

                cc = get_kbrd(1)
                
             endif

             mkhdr, header, objflux, /extend
             sxdelpar, header, 'COMMENT'

             sxaddpar, header, 'RA', zmatch.ra
             sxaddpar, header, 'DEC', zmatch.dec
             sxaddpar, header, 'GALAXY', galaxy
             sxaddpar, header, 'Z', float(zmatch.z_obj), ' redshift'
             sxaddpar, header, 'Z_ERR', float(0.0), ' redshift error'
             sxaddpar, header, 'CRVAL1', minwave, ' wavelength at CRPIX1'
             sxaddpar, header, 'CRPIX1', 1.0, ' reference pixel number'
             sxaddpar, header, 'CDELT1', dwave, ' dispersion [Angstrom/pixel]'
             sxaddpar, header, 'CTYPE1', 'LINEAR', ' projection type'
             sxaddpar, header, 'CD1_1', dwave, ' dispersion [Angstrom/pixel]'

             if keyword_set(wfits) and (zmatch.z_obj gt 0.0) and (zmatch.z_obj lt 1.0) then begin

                outfile = galaxy+'.fits'

                splog, 'Writing '+outpath+outfile+'.'
                mwrfits, objflux, outpath+outfile, header, /create
                mwrfits, objferr, outpath+outfile
                spawn, ['gzip -f '+outpath+outfile], /sh
                
             endif

          endif
          
       endfor 

    endfor 

stop
    
return
end
