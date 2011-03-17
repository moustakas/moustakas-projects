pro sings_write_extranucs, final, airmass=airmass, catalog=catalog, visualize=visualize, offset=offset
; jm06mar24uofa

    all = rsex('sings_extranucs.txt')
    moire = rsex('sings_prescott_extranucs.txt')
    nall = n_elements(all)

    final = struct_addtags(all,replicate({fluxha: -999.0, flux24: -999.0, regionx: -999, $
      regiony: -999, vmag: 99, object: '', priority: 1L, id: ''},nall))

    for i = 0L, nall-1L do begin
       match = where((all[i].ra eq moire.ra) and (all[i].dec eq moire.dec),nmatch)
       if (nmatch eq 1L) then begin
          final[i].hii = 'P'+string(moire[match].id,format='(I0)')
          final[i].fluxha = moire[match].fluxha;*1E-3*1D-23*2.99D18/6563.0^2*1D16
          final[i].flux24 = moire[match].flux24
          final[i].regionx = moire[match].regionx
          final[i].regiony = moire[match].regiony
       endif
    endfor

    final = final[sort(final.ra)]
    keep = where((im_hms2dec(final.ra) gt 5.0) and (im_hms2dec(final.ra) lt 24.0),nfinal)
    final = final[keep]
;   final.id = final.object
    final.id = string(lindgen(nfinal)+1L,format='(I2.2)')
    final.object = final.id+'/'+repstr(strtrim(final.galaxy,2),'NGC','N')+'/'+strtrim(final.hii,2)+'/'+strtrim(final.selection,2)

    mwrfits, final, 'extranucs.fits', /create
    
    galaxy = strtrim(final.galaxy,2)
    unique = uniq(galaxy)
    ugalaxy = galaxy[unique]
    ngalaxy = n_elements(ugalaxy)

    ufinal = final[unique]

    if keyword_set(airmass) then begin
       airmass_plots, '2006-05-21', final.ra, final.dec, object=final.object, /bigpost, psname='extranuc_06may.ps'
       airmass_plots, '2006-05-21', final[unique].ra, final[unique].dec, object=final[unique].galaxy, $
         /bigpost, psname='extranuc_galaxy_06mar.ps'
    endif

    if keyword_set(catalog) then begin

       openw, unit, 'extranuc_06mar_catalog.txt', /get_lun
       printf, unit, '# John Moustakas -- 2006 Mar 25 -- jmoustakas@as.arizona.edu'
       printf, unit, '# ID              ', 'RA', 'DEC', 'RA_PM', 'DEC_PM', 'Vmag', $
         'P', 'Equinox', format='(A15,A13,A13,A8,A8,A6,A3,A9)' 
       printf, unit, '#---------------------------------------------------------' + $
         '-----------------'

       for i = 0L, nfinal-1L do begin
          printf, unit, final[i].id + '               ', final[i].ra, final[i].dec, $
            '+00.000', '+00.000', final[i].vmag, final[i].priority, 'J2000.0', $
            format='(A15,A13,A13,A8,A8,F6.2,I3,A9)'
       endfor
       free_lun, unit

    endif

    if keyword_set(visualize) then begin

       im_window, 0, xratio=0.6, /square
       
       rad = 8.0/60.0 ; [arcsec]
       
       dsspath = sings_path()+'DSS/'
       dssfits = dsspath+strlowcase(ugalaxy)+'.fits.gz'

       for i = 0L, ngalaxy-1L do begin

          match = where(galaxy eq ugalaxy[i],nmatch)
          
          dssimage = readfits(dssfits[i],hdss,/silent)
          gsssextast, hdss, astr
          
          if (file_test(dssfits[i],/regular) eq 0L) then begin
             splog, 'DSS image for '+ugalaxy[i]+' not found!'
             return
          endif

          imsize = size(dssimage,/dimension)
          xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

          xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
          ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

          xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
          yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

          img = logscl(dssimage,exponent=1.5,negative=keyword_set(postscript),omin=35,omax=255)

; hii regions          

          gsssadxy, astr, 15.0*im_hms2dec(final[match].ra), im_hms2dec(final[match].dec), xhii, yhii
          xhii = (xhii - xcen)*xpixscale & yhii = (yhii - ycen)*ypixscale

;         imgxrange = max(abs(xhii))*[-1.2,1.2] & imgyrange = max(abs(yhii))*[-1.2,1.2]
          imgyrange = minmax(yaxis) & imgxrange = minmax(xaxis)
          
; plot the image          
          
          plotimage, img, /preserve_aspect, position=pos, /normal, imgxrange=imgxrange, $
            imgyrange=imgyrange, charsize=1.8, charthick=postthick, xthick=postthick, $
            ythick=postthick, xtitle=textoidl('\Delta\alpha [arcmin]'), $
            ytitle=textoidl('\Delta\delta [arcmin]'), title=name
          legend, ugalaxy[i], /left, /top, box=0, charsize=1.5, charthick=2.0

; overplot the HII regions       

          for ii = 0L, nmatch-1L do begin
             tvcircle, rad, xhii[ii], yhii[ii], color=djs_icolor('red'), $
               thick=2.0, /data, line=0
          endfor

          cc = get_kbrd(1)
          
       endfor

    endif

    if keyword_set(offset) then begin

       offset = replicate({galaxy: '', star_ra: '', star_dec: '', rmag: 99.0, bmag: 99.0},ngalaxy)
       offset.galaxy = ufinal.galaxy
;      for i = 13L, ngalaxy-1L do begin
       for i = 0L, ngalaxy-1L do begin
          ra = im_hms2dec(ufinal[i].ra)*15.0 & dec = im_hms2dec(ufinal[i].dec)
          info = queryusno([ra,dec],20.0,magrange=[5,15])
;         info = queryvizier('USNO',[ra,dec],15.0)
          if (size(info,/type) eq 8L) then begin
             good = where((info.b_mag gt 5.0) and (info.b_mag lt 13.0),ngood)
             niceprint, info[good].id, im_dec2hms(info[good].ra/15.0,/col), im_dec2hms(info[good].dec,/col), info[good].b_mag
;            ngood = n_elements(info) & good = lindgen(ngood)
;            good = where(info.recno ne 0L,ngood)
             if (ngood ne 0L) and (tag_exist(info,'RA')) then begin
                star_ra = info[good].ra & star_dec = info[good].dec
;               star_ra = im_hms2dec(info[good].ra)*15.0 & star_dec = im_hms2dec(info[good].dec)
                diff = djs_diff_angle(ra,dec,star_ra,star_dec,units='degrees')
                mindiff = min(diff,indx)
                offset[i].star_ra = im_dec2hms(info[good[indx]].ra/15.0,/colon)
                offset[i].star_dec = im_dec2hms(info[good[indx]].dec,/colon)
                offset[i].rmag = info[good[indx]].r_mag
                offset[i].bmag = info[good[indx]].b_mag
             endif else splog, 'No offset star found for '+ufinal[i].galaxy
          endif
          splog, strtrim(offset[i].galaxy,2)+': press any key to continue.'
          cc = get_kbrd(1)
       endfor

       struct_print, offset, /no_head
       mwrfits, offset, 'extranucs_offset_stars.fits', /create
       struct_print, offset, file='extranucs_offset_stars.txt'

stop       
       
    endif
    
return
end
    
