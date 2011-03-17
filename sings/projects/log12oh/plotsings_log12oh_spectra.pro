pro plotsings_log12oh_spectra, all=all, postscript=postscript
; jm10mar11ucsd - generate the visualizations of the galaxies and the
; spectra 

    if keyword_set(all) then postscript = 1
    
; grep "LL/SL nuc" mips_irac_irs_schedule_05mar22.txt > irs_schedule_05mar22.txt
; mpage -2 -l -W140 irs_schedule_05mar22.txt > irs_schedule_05mar22.ps

; Y is "along" the slit, X is perpendicular to the slit. 

; slit_length = 3.5*60; slit_overlap = 40.0; the first slit offset is
; first_offset = (slit_length-slit_overlap)/2.0 = 85"; subsequently,
; slit_offsets = first_offset + (lindgen(noffsets)+1)*slit_length - slit_overlap/2.0

    obspath = sings_path(/observing)
    dsspath = sings_path(/dss)
    analysis_path = sings_path(/analysis)
    spec1dpath = sings_path(/spec1d)
    paperpath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'

    sings = sings_read_info()

    galaxy = strtrim(sings.galaxy,2)
    nicegalaxy = strtrim(sings.nice_galaxy,2)
    ngalaxy = n_elements(galaxy)

; read the PPXF fitting results
    allnuclear = read_sings_gandalf(/nuclear)
    alldrift20 = read_sings_gandalf(/drift20)
    alldrift56 = read_sings_gandalf(/drift56)

    allnuclear_specfit = read_sings_gandalf_specfit(allnuclear,/nuclear)
    alldrift20_specfit = read_sings_gandalf_specfit(alldrift20,/drift20)
    alldrift56_specfit = read_sings_gandalf_specfit(alldrift56,/drift56)
    
; read the April spectroscopy sample, driftscan parameters, and the
; NED data; SAMPLE and DATA both contain the NED galaxy centers; the
; IRS centers, when available (see IRS_GETDSS), are in the DSS
; headers; also read the IRS parameters for other quantities of
; interest 
    irs = rsex(obspath+'irs_ll_params.txt')
    dssfits = dsspath+strcompress(strlowcase(nicegalaxy),/remove)+'.fits.gz'
    
; FOR THE PAPER!
    
    if keyword_set(all) then begin
       
       nperpage = 6L
       npage = ceil(ngalaxy/float(nperpage))

       psname = 'FIG_LOG12OH/'+string(lindgen(ngalaxy)+1L,form='(I2.2)')+'.'+strlowcase(galaxy)+'.eps'
       suffix = strarr(ngalaxy) & suffix[lindgen(ngalaxy/(nperpage/2L))*(nperpage/2L)+2] = ' \\'

;      npage = 3L               ; TEMPORARY!

       label = string(lindgen(npage)+1,format='(I2.2)')
       label[0] = 'first'
       label[npage-1L] = 'last'

       openw, lun, sings_path(/papers)+'log12oh/spectra.txt', /get_lun
       for ipage = 0L, npage-1L do begin
          x1 = nperpage*ipage
          x2 = (nperpage*(ipage+1))<(ngalaxy-1L)
          these = [reverse(psname[x1:x1+nperpage/2L-1L]),reverse(psname[nperpage/2L:x2-1L])]
          suffix = ['','',' \\','','','']
          printf, lun, '\begin{figure}[!h]'
          printf, lun, '\begin{center}'
          niceprintf, lun, '\includegraphics[scale=0.42,angle=0]{'+these+'}'+suffix
          printf, lun, '\caption{Caption.\label{fig:spectra_'+label[ipage]+'}}'
          printf, lun, '\end{center}'
          printf, lun, '\end{figure}'
          printf, lun, ' '
          printf, lun, '\clearpage'
       endfor
       free_lun, lun

    endif

; change these values with care; they were picked by looking at the
; output from PLOTIMAGE after setting PRESERVE_ASPECT=1; set the size
; of the DSS image to 0.41875 in "normal" units
    xpage = 11.0

    xmargin1 = 0.85
    xmargin2 = 0.15
    ymargin1 = 0.6
    ymargin2 = 0.1
    xdss = 4.2
    ydss = xdss
    xspec = xpage-xdss-xmargin1-xmargin2
    yspec = ydss/3.0
    
    ypage = ymargin1 + ymargin2 + ydss

    pos = fltarr(4,4)
    pos[*,0] = [xmargin1,ymargin1+2*yspec,xmargin1+xspec,ymargin1+3*yspec]
    pos[*,1] = [xmargin1,ymargin1+1*yspec,xmargin1+xspec,ymargin1+2*yspec]
    pos[*,2] = [xmargin1,ymargin1+0*yspec,xmargin1+xspec,ymargin1+1*yspec]
    pos[*,3] = [xpage-xdss-xmargin2,ymargin1,xpage-xmargin2,ymargin1+ydss]
    
    pos[[0,2],*] = pos[[0,2],*] / xpage
    pos[[1,3],*] = pos[[1,3],*] / ypage
    
    if keyword_set(postscript) then begin
       im_plotfaves, /postscript
       if keyword_set(all) then begin
;         dfpsplot, paperpath+'sings_log12oh_spectra.ps', /landscape, xsize=xpage, $
;           ysize=ypage, xoffset=((8.5-ypage)/2.0)>0.0, yoffset=(11.0-(11.0-xpage)/2.0)>0.0, $
;           /color
          set_plot, 'PS'
          device, file=paperpath+'sings_log12oh_spectra.ps', /landscape, xsize=xpage, ysize=ypage, $
            xoffset=((8.5-ypage)/2.0)>0.0, yoffset=(11.0-(11.0-xpage)/2.0)>0.0, /inch, /color, bits=24
       endif
    endif else begin
       im_window, 0, xratio=0.7, yratio=0.7
    endelse

    charsize1 = 1.2
    
;   for i = 0, 1 do begin
;   for i = 15, ngalaxy-1L do begin
    for i = 0, ngalaxy-1 do begin

       if keyword_set(postscript) then begin
          if (not keyword_set(all)) then begin
             prefix = string(i+1,format='(I2.2)')+'.'
             set_plot, 'PS'
             device, file=paperpath+prefix+strlowcase(galaxy[i])+'.eps', /landscape, xsize=xpage, ysize=ypage, $
               xoffset=((8.5-ypage)/2.0)>0.0, yoffset=(11.0-(11.0-xpage)/2.0)>0.0, /inch, $
               /color, /encapsulated, bits=24
          endif
       endif

;      zobj = sings[i].z
       
; ###########################################################################       
; visualization the DSS image 
; ###########################################################################       
       
       dssimage = readfits(dssfits[i],hdss,/silent)
       gsssextast, hdss, astr

       irs_ra = repstr(strtrim(sxpar(hdss,'OBJCTRA'),2),' ',':')
       irs_dec = repstr(strtrim(sxpar(hdss,'OBJCTDEC'),2),' ',':')
       
       imsize = size(dssimage,/dimension)
       xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0

       xpixscale = (astr.pltscl*1E-3*astr.xsz)/60 ; [arcmin/pixel]
       ypixscale = (astr.pltscl*1E-3*astr.ysz)/60 ; [arcmin/pixel]

       xaxis = (findgen(xsize)-xcen)*xpixscale ; centered on the image [arcsec]
       yaxis = (findgen(ysize)-ycen)*ypixscale ; centered on the image [arcsec]

       img = logscl(dssimage,exponent=1.5,negative=keyword_set(postscript),omin=35,omax=255)
       plotimage, img, preserve_aspect=1, position=pos[*,3], /normal, imgxrange=minmax(xaxis), $
         imgyrange=minmax(yaxis), charsize=1.8, xtitle='', ytitle='', /noaxes
       legend, textoidl(nicegalaxy[i]), /right, /top, box=0, charsize=1.5, $
         textcolor=djs_icolor('black'), /normal, margin=0

       if (1.2*sings[i].d25_maj le 3.0) then begin
          xbar = 30.0
          barlabel = '30 arcsec'
       endif
       if (1.2*sings[i].d25_maj gt 3.0) and (1.2*sings[i].d25_maj lt 15.0) then begin
          xbar = 60.0
          barlabel = '1 arcmin'
       endif
       if (1.2*sings[i].d25_maj ge 15.0) then begin
          xbar = 180.0
          barlabel = '5 arcmin'
       endif

       im_oplot_box, 0.0, xbar/60.0, 90.0, xoffset=-xsize*xpixscale*0.8/2.0, $
         yoffset=+ysize*ypixscale*0.7/2.0, line=0, $
         corners=corners, /noplot
       oplot, [corners[0,1],corners[0,0]], corners[1,0]*[1,1], $
         line=0, thick=!p.thick+2
       xyouts, (corners[0,1]-corners[0,0])/2.0+corners[0,0], corners[1,0]*1.15, $
         barlabel, charsize=1.0, align=0.5, /data

; overplot the RC3 ellipse
       gsssadxy, astr, 15.0*im_hms2dec(sings[i].ra), im_hms2dec(sings[i].dec), xrc3, yrc3
       xrc3 = (xrc3 - xcen)*xpixscale & yrc3 = (yrc3 - ycen)*ypixscale
       
       if (sings[i].posangle gt -900.0) then tvellipse, sings[i].d25_maj/2.0, $
         sings[i].d25_min/2.0, xrc3, yrc3, 90+sings[i].posangle, color=djs_icolor('purple'), $
         /data, line=2

; check if there is an IRS spectrum
       match, strlowcase(strcompress(nicegalaxy[i],/remove)), $
         strlowcase(strtrim(irs.nice_galaxy,2)), indx1, indx2
       if (indx1[0] ne -1L) then begin
          irsinfo = irs[indx2]
          irs_ra = irsinfo.ra & irs_dec = irsinfo.dec
       endif else begin
          irsinfo = -1L
          irs_ra = sings[i].ra & irs_dec = sings[i].dec
       endelse

; visualize the spectroscopic apertures

       match, strlowcase(galaxy[i]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
       if (indx1[0] eq -1L) then $
         specinfo = {drift56: 0L, drift20: 0L, nuclear: 0L} else $
           specinfo = sings[indx2]

       if specinfo.drift56 then begin

          if strmatch(galaxy[i],'*NGC3034*',/fold) then $
            gsssadxy, astr, 15.0*im_hms2dec(irsinfo.ra), im_hms2dec(irsinfo.dec), x56, y56 else $ ; NOTE coordinates
              gsssadxy, astr, 15.0*im_hms2dec(sings[i].ra), im_hms2dec(sings[i].dec), x56, y56
          x56 = (x56 - xcen)*xpixscale
          y56 = (y56 - ycen)*ypixscale

          im_oplot_box, specinfo.drift56_scan/60.0, specinfo.drift56_ap/60.0, $
            specinfo.drift56_posangle, /noplot, corners=corners

          indx = [0,1,2,3,0]
          djs_oplot, [corners[0,indx],corners[0,indx+1]]+x56, $
            [corners[1,indx],corners[1,indx+1]]+y56, line=0, $
            color=djs_icolor('red')
          
       endif

       if specinfo.drift20 then begin

          gsssadxy, astr, 15.0*im_hms2dec(sings[i].ra), im_hms2dec(sings[i].dec), x20, y20
          x20 = (x20 - xcen)*xpixscale
          y20 = (y20 - ycen)*ypixscale

          im_oplot_box, specinfo.drift20_scan/60.0, specinfo.drift20_ap/60.0, $
            specinfo.drift20_posangle, /noplot, corners=corners

          indx = [0,1,2,3,0]
          djs_oplot, [corners[0,indx],corners[0,indx+1]]+x20, $
            [corners[1,indx],corners[1,indx+1]]+y20, line=0, $
            color=djs_icolor('blue')

       endif
          
; overlay the IRS aperture

       if (size(irsinfo,/type) eq 8L) then begin

          gsssadxy, astr, 15.0*im_hms2dec(irsinfo.ra), im_hms2dec(irsinfo.dec), xirs, yirs
          xirs = (xirs - xcen)*xpixscale
          yirs = (yirs - ycen)*ypixscale

          im_oplot_box, irsinfo.width, irsinfo.length_2x, im_angle_format(irsinfo.pa), $
            /noplot, corners=corners

          indx = [0,1,2,3,0]
          djs_oplot, [corners[0,indx],corners[0,indx+1]]+xirs, $
            [corners[1,indx],corners[1,indx+1]]+yirs, line=1, $
            color=djs_icolor('yellow')

       endif

; ###########################################################################       
; now visualization the spectra
; ###########################################################################       

       factor = 1.0
       sigrej = 4.0
       scale = 1E15
       ytitle = 'Flux (10^{-15} '+flam_units()+')'
       xtitle = 'Rest Wavelength (\AA)'

       match, strlowcase(galaxy[i]), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2
       if (indx1[0] eq -1L) then $
         specinfo = {drift56: 0L, drift20: 0L, nuclear: 0L} else $
           specinfo = sings[indx2]

       if specinfo.drift56 then begin

          this = where(strtrim(alldrift56.galaxy,2) eq galaxy[i])
          good = where(alldrift56_specfit[this].wave gt 0.0)
          wave = exp(alldrift56_specfit[this].wave[good])
          flux = scale*alldrift56_specfit[this].flux[good]
          linefit = scale*alldrift56_specfit[this].linefit[good]
          continuum = scale*alldrift56_specfit[this].continuum[good]
          
;         specdata = rd1dspec(strtrim(specinfo.drift56_file,2),datapath=spec1dpath,/silent)
;         wave = specdata.wave
;         flux = scale*specdata.spec
          refwave = wave
;         emask = emission_mask(wave*(1+zobj),z=zobj)
          
          case strcompress(strlowcase(galaxy[i]),/remove) of
             'ngc0925': yscale = [0.95,1.1]
             'ngc1097': yscale = [0.7,0.9]
             'ngc1705': yscale = [1.0,1.3]
             'ngc2403': yscale = [0.95,1.15]
             'm81dwa':  yscale = [0.9,1.2]
             'ddo053':  yscale = [0.9,1.2]
             'ngc2915': yscale = [0.95,1.25]
             'ngc2976': yscale = [0.95,1.15]
             else: yscale = [0.95,1.05]
          endcase

          stats = im_stats((continuum+linefit)<flux,sigrej=sigrej)
;         stats = im_stats(continuum+linefit,sigrej=sigrej)
;         stats = im_stats(flux[where(emask)],sigrej=sigrej)
          yrange = [yscale[0]*stats.minrej,stats.maxrej*yscale[1]]
          xrange = minmax(wave)
          
          djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xtickname=replicate(' ',10), $
            xtitle='', ytitle='', charsize=charsize1, $
            yrange=yrange, xrange=xrange, position=pos[*,0], yminor=3
          djs_oplot, wave, flux, ps=10, thick=1.0, color='grey' ;, color='red'
          legend, 'Radial Strip', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data

;         spec = read_sings_specfit(galaxy[i],/drift56,/silent)
;         djs_oplot, spec[*,0]*(1+zobj), factor*scale*spec[*,2]/(1+zobj), thick=postthick2
          djs_oplot, wave, continuum, psym=10, thick=2.0

       endif else begin

          djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xtickname=replicate(' ',10), $
            xtitle='', ytitle='', charsize=charsize1, $
            position=pos[*,0], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
;         xyouts, (pos[2,0]-pos[0,0])/2.0+pos[0,0], (pos[3,0]-pos[1,0])/2.0+pos[1,0], $
;           'Spectrum Unavailable', align=0.5, charthick=postthick, /normal
          legend, 'Radial Strip Spectrum Unavailable', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data
          
       endelse

       if specinfo.drift20 then begin

          this = where(strtrim(alldrift20.galaxy,2) eq galaxy[i])
          good = where(alldrift20_specfit[this].wave gt 0.0)
          wave = exp(alldrift20_specfit[this].wave[good])
          flux = scale*alldrift20_specfit[this].flux[good]
          linefit = scale*alldrift20_specfit[this].linefit[good]
          continuum = scale*alldrift20_specfit[this].continuum[good]

;         specdata = rd1dspec(strtrim(specinfo.drift20_file,2),datapath=spec1dpath,/silent)
;         wave = specdata.wave
;         flux = scale*specdata.spec

          if (specinfo.drift56 eq 0) then refwave = wave else begin
             flux = interpol(flux,wave,refwave)
             linefit = interpol(linefit,wave,refwave)
             continuum = interpol(continuum,wave,refwave)
             wave = refwave
          endelse
          
;         emask = emission_mask(wave*(1+zobj),z=zobj)

          case strcompress(strlowcase(galaxy[i]),/remove) of
             'ngc0925':   yscale = [0.95,1.15]
             'ngc1097':   yscale = [0.7,0.9]
             'ngc1705':   yscale = [0.95,1.1]
             'm81dwa':    yscale = [0.9,1.2]
             'ddo053':    yscale = [0.9,1.2]
             'holmbergi': yscale = [0.9,1.25]
             'ngc2915':   yscale = [0.95,1.25]
             'ngc2976':   yscale = [0.95,1.15]
             'ngc3049':   yscale = [0.95,1.15]
             else: yscale = [0.95,1.05]
          endcase

          stats = im_stats((continuum+linefit)<flux,sigrej=sigrej)
;         stats = im_stats(continuum+linefit,sigrej=sigrej)
;         stats = im_stats(flux[where(emask)],sigrej=sigrej)
          yrange = [yscale[0]*stats.minrej,stats.maxrej*yscale[1]]
          xrange = minmax(wave)
          
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xtickname=replicate(' ',10), $
            xtitle='', ytitle=textoidl(ytitle)+'!c', charsize=charsize1, /noerase, $
            yrange=yrange, xrange=xrange, position=pos[*,1], yminor=3
          djs_oplot, wave, flux, ps=10, thick=1.0, color='grey'; color='blue'

          legend, 'Circumnuclear', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data

;         spec = read_sings_specfit(galaxy[i],/drift20,/silent)
;         djs_oplot, spec[*,0]*(1+zobj), factor*scale*spec[*,2]/(1+zobj), thick=postthick2
          djs_oplot, wave, continuum, psym=10, thick=2.0

       endif else begin

          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xtickname=replicate(' ',10), $
            xtitle='', ytitle='', charsize=charsize1, /noerase, $
            position=pos[*,1], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
;         xyouts, (pos[2,1]-pos[0,1])/2.0+pos[0,1], (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
;           'Spectrum Unavailable', align=0.5, charthick=postthick, /normal
          legend, 'Circumnuclear Spectrum Unavailable', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data
          
       endelse
       
       if specinfo.nuclear then begin

          this = where(strtrim(allnuclear.galaxy,2) eq galaxy[i])
          good = where(allnuclear_specfit[this].wave gt 0.0)
          wave = exp(allnuclear_specfit[this].wave[good])
          flux = scale*allnuclear_specfit[this].flux[good]
          linefit = scale*allnuclear_specfit[this].linefit[good]
          continuum = scale*allnuclear_specfit[this].continuum[good]

;         specdata = rd1dspec(strtrim(specinfo.nuclear_file,2),datapath=spec1dpath,/silent)
;         wave = specdata.wave
;         flux = scale*specdata.spec

          if (specinfo.drift56 eq 0) and (specinfo.drift20 eq 0) then refwave = wave else begin
             flux = interpol(flux,wave,refwave)
             linefit = interpol(linefit,wave,refwave)
             continuum = interpol(continuum,wave,refwave)
             wave = refwave
          endelse
          
;         emask = emission_mask(wave*(1+zobj),z=zobj)

          case strcompress(strlowcase(galaxy[i]),/remove) of
             'ngc0024': yscale = [0.95,1.1]
             'ngc0925': yscale = [0.95,1.15]
             'ngc1097': yscale = [0.6,0.9]
             'ngc1705': yscale = [1.0,1.1]
             'ngc2915': yscale = [0.95,1.25]
             'ngc2976': yscale = [0.95,1.15] 
            else: yscale = [0.95,1.05]
          endcase

          stats = im_stats((continuum+linefit)<flux,sigrej=sigrej)
;         stats = im_stats(continuum+linefit,sigrej=sigrej)
;         stats = im_stats(flux[where(emask)],sigrej=sigrej)
          yrange = [yscale[0]*stats.minrej,stats.maxrej*yscale[1]]
          xrange = minmax(wave)
          
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xtitle=xtitle, $
            ytitle='', charsize=charsize1, /noerase, yrange=yrange, xrange=xrange, $
            position=pos[*,2], yminor=3
          djs_oplot, wave, flux, ps=10, thick=1.0, color='grey'; color='dark green'

          legend, 'Nuclear', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data
          
;         spec = read_sings_specfit(galaxy[i],/nuclear,/silent)
;         djs_oplot, reform(spec[*,0])*(1+zobj), factor*scale*reform(spec[*,2])/(1+zobj), thick=postthick2
          djs_oplot, wave, continuum, psym=10, thick=2.0

       endif else begin

          if (specinfo.drift56 eq 0) and (specinfo.drift20 eq 0) then begin
             
             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xtickname=replicate(' ',10), $
               xtitle='', ytitle='', charsize=charsize1, /noerase, $
               position=pos[*,2], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
             
          endif else begin

             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, $
               xtitle=xtitle, ytitle='', charsize=charsize1, /noerase, $
               position=pos[*,2], ytickname=replicate(' ',10), yticklen=1E-10, xticklen=1E-10
             
          endelse

;         xyouts, (pos[2,2]-pos[0,2])/2.0+pos[0,2], (pos[3,2]-pos[1,2])/2.0+pos[1,2], $
;           'Spectrum Unavailable', align=0.5, charthick=postthick, /normal
          legend, 'Nuclear Spectrum Unavailable', /left, /top, box=0, charsize=charsize1, $
            position=[!x.crange[0]*1.07,!y.crange[1]*0.95], /data

       endelse

;      if (not keyword_set(postscript)) then cc = get_kbrd(1)
       if keyword_set(postscript) then begin
          if (keyword_set(all) eq 0) then dfpsclose
       endif else cc = get_kbrd(1)

    endfor 

    if keyword_set(all) then begin
       dfpsclose
       im_plotfaves
    endif
    
stop       

return
end    
