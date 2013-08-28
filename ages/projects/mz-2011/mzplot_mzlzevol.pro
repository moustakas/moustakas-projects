pro mzplot_mzlzevol, ps=ps
; jm09mar27nyu - plot the evolution of the LZ and MZ relations
; jm10oct10ucsd - major update    

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; read the data    
    zbins = mz_zbins(nzbins)
    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)

    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)
    agesohdust = read_mz_sample(/mzhii_log12oh)
    agesohnodust = read_mz_sample(/mzhii_log12oh,/nodust)

    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1)

; --------------------------------------------------
; Figure 18 - AGES - LZ evolution
    xrange = [-17.2,-23.2]
;   yrange = [8.35,9.37] ; [8.3,9.39]
    yrange = [8.4,9.3] ; [8.3,9.39]

    localline = 0
    localcolor = 'black'
    evolline = 5
    levels = [0.25,0.5,0.75,0.9]

    calibcolor = ['forest green','dodger blue','firebrick']
    below_calibpsym = [5,6,9]
    above_calibpsym = [17,15,16]
    calibsymsize = [1.5,1.8,1.5]*0.9
    
    for ii = 0, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
             ohrange1 = [8.3,9.35]
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
             ohrange1 = [8.3,9.2]
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
             ohrange1 = [8.55,9.35]
          end
       endcase
       ytitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
       xtitle1 = mzplot_mbtitle()
       
       lzlocal = mrdfits(mzpath+'lzlocal_sdss_ews_'+calib+'.fits.gz',1)
       lzevol = mrdfits(mzpath+'lzevol_'+calib+'.fits.gz',1)
       mzevol = mrdfits(mzpath+'mzevol_'+calib+'.fits.gz',1)
       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit,/errcut) ; apply an error cut for the contours
       
       psfile = pspath+'lzevol_'+calib+suffix
       im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
         height=2.7*[1,1,1]

       zbins = mz_zbins(nzbins)
       for iz = 0, nzbins-1 do begin
          these = where((ainfo.z ge zbins[iz].zlo) and $
            (ainfo.z lt zbins[iz].zup),nthese)
          xx = ainfo.mb_ab[these]
          yy = ainfo.oh[these]
          ww = ainfo.weight[these]
          ee = ainfo.oh_err[these]
          
          if odd(iz) then begin
             ytitle = '' & ytickname = replicate(' ',10)
          endif else begin
             ytitle = ytitle1 & delvarx, ytickname
          endelse
          if (iz lt 4) then begin
             xtitle = '' & xtickname = replicate(' ',10)
          endif else delvarx, xtickname 

          mzplot_scatterplot, xx, yy, weight=ww, noerase=(iz gt 0), $
            position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
            xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
            ytickname=ytickname, _extra=extra, ccolor=im_color('grey30',201), $
            /nogrey, levels=levels, xtickinterval=2
          djs_oplot, limits.mblimit_50[iz]*[1,1], yrange+[-0.1,+0.1], line=1, thick=8

;;          above = where(mzevol.complete_bymass[*,iz] eq 1 and mzevol.ohmean_bymass[*,iz] gt -900,nabove)
;;          if (nabove ne 0) then begin
;;             oploterror, mzevol.medmass_bymass[above,iz], mzevol.ohmean_bymass[above,iz], $
;;               mzevol.ohmean_bymass_err[above,iz], psym=symcat(above_calibpsym[ii],thick=5), $
;;               symsize=calibsymsize[ii], color=im_color(calibcolor[ii]), $
;;               errcolor=im_color(calibcolor[ii])
;;          endif
;;          
;;          below = where(mzevol.complete_bymass[*,iz] eq 0 and mzevol.ohmean_bymass[*,iz] gt -900,nbelow)
;;          if (nbelow ne 0) then begin
;;             oploterror, mzevol.medmass_bymass[below,iz], mzevol.ohmean_bymass[below,iz], $
;;               mzevol.ohmean_bymass_err[below,iz], psym=symcat(below_calibpsym[ii],thick=5), $
;;               symsize=calibsymsize[ii], color=im_color(calibcolor[ii]), $
;;               errcolor=im_color(calibcolor[ii])
;;          endif
          
          legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
            '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
            /left, /top, box=0, margin=0, charsize=1.4

; local LZ          
;         oplot_lzfit, lzlocal.coeff, z=zbins[0].zbin, qz0=lzevol.qz0, $
;           linestyle=localline, linecolor=localcolor

          oplot_lzfit, lzevol.coeff, linestyle=localline, $
            linecolor=localcolor, z=zbins[0].zbin, qz0=lzevol.qz0, /evol
          if iz gt 0 then begin
             oplot_lzfit, lzevol.coeff, linestyle=evolline, $
               linecolor=calibcolor[ii], z=zbins[iz].zbin, qz0=lzevol.qz0, /evol
             coeff = lzevol.coeff
             coeff[2] = mzavg.mzevol_coeffs_r0zero[4,ii]
             oplot_lzfit, coeff, linestyle=3, $
               linecolor='midnight blue', z=zbins[iz].zbin, qz0=lzevol.qz0, /evol
          endif
          xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, xtitle1
       endfor 
       im_plotconfig, /psclose, psfile=psfile
    endfor

; print, (abs(mzavg.p0r0zero_avg)-abs(mzavg.lz_s0_avg))/mzavg.lz_slope_avg

; --------------------------------------------------
; Figure 13 - AGES - MZ evolution
    xrange = [8.8,11.3]
    yrange = [8.4,9.3] ; [8.3,9.39]
;   xrange = [8.4,11.7]
;   yrange = [8.35,9.37] ; [8.3,9.39]
    maxis = range(xrange[0]+0.1,xrange[1]-0.1,100)

    localline = 0
    localcolor = 'black'
    evolline = 5
    evolcolor = 'firebrick'
    levels = [0.25,0.5,0.75,0.9]

    calibcolor = ['forest green','dodger blue','firebrick']
    below_calibpsym = [5,6,9]
    above_calibpsym = [17,15,16]
    calibsymsize = [1.5,1.8,1.5]*0.9
    
    for ii = 0, 2 do begin
       case ii of
          0: begin
             t04 = 1 & m91 = 0 & kk04 = 0
             calib = 't04'
             ohrange1 = [8.3,9.35]
          end
          1: begin
             t04 = 0 & m91 = 1 & kk04 = 0
             calib = 'm91'
             ohrange1 = [8.3,9.2]
          end
          2: begin
             t04 = 0 & m91 = 0 & kk04 = 1
             calib = 'kk04'
             ohrange1 = [8.55,9.35]
          end
       endcase
       ytitle1 = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
       xtitle1 = mzplot_masstitle()
       
       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
       mzevol = mrdfits(mzpath+'mzevol_'+calib+'.fits.gz',1)

;      ainfo = mzlz_grab_info(agesohnodust,agesancillary,agesmass,$
;        t04=t04,m91=m91,kk04=kk04,/nolimit,/flux,zmin=0.05,zmax=0.15)
       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
         t04=t04,m91=m91,kk04=kk04,/nolimit,/errcut) ; apply an error cut for the contours
       
       psfile = pspath+'mzevol_'+calib+suffix
       im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
         height=2.7*[1,1,1]

       zbins = mz_zbins(nzbins)
       for iz = 0, nzbins-1 do begin
          these = where((ainfo.z ge zbins[iz].zlo) and $
            (ainfo.z lt zbins[iz].zup),nthese)
          xx = ainfo.mass[these]
          yy = ainfo.oh[these]
          ww = ainfo.weight[these]
          ee = ainfo.oh_err[these]
          
          if odd(iz) then begin
             ytitle = '' & ytickname = replicate(' ',10)
          endif else begin
             ytitle = ytitle1 & delvarx, ytickname
          endelse
          if (iz lt 4) then begin
             xtitle = '' & xtickname = replicate(' ',10)
          endif else delvarx, xtickname 

          mzplot_scatterplot, xx, yy, weight=ww, noerase=(iz gt 0), $
            position=pos[*,iz], xrange=xrange, yrange=yrange, xsty=1, ysty=1, $
            xtitle=xtitle, ytitle=ytitle, xtickname=xtickname, $
            ytickname=ytickname, _extra=extra, ccolor=im_color('grey30',201), $
            /nogrey, levels=levels, xtickinterval=1
;         djs_oplot, limits.masslimit_50[iz]*[1,1], yrange+[+0.1,-0.1], line=1, thick=8

          above = where(mzevol.complete_bymass[*,iz] eq 1 and mzevol.ohmean_bymass[*,iz] gt -900,nabove)
          if (nabove ne 0) then begin
             oploterror, mzevol.medmass_bymass[above,iz], mzevol.ohmean_bymass[above,iz], $
               mzevol.ohmean_bymass_err[above,iz], psym=symcat(above_calibpsym[ii],thick=5), $
               symsize=calibsymsize[ii], color=im_color(calibcolor[ii]), $
               errcolor=im_color(calibcolor[ii])
          endif
          
          below = where(mzevol.complete_bymass[*,iz] eq 0 and mzevol.ohmean_bymass[*,iz] gt -900,nbelow)
          if (nbelow ne 0) then begin
             oploterror, mzevol.medmass_bymass[below,iz], mzevol.ohmean_bymass[below,iz], $
               mzevol.ohmean_bymass_err[below,iz], psym=symcat(below_calibpsym[ii],thick=5), $
               symsize=calibsymsize[ii], color=im_color(calibcolor[ii]), $
               errcolor=im_color(calibcolor[ii])
          endif
          
;         abin = im_medxbin(xx,yy,0.15,weight=ww/ee^2,/verbose,minpts=3)
;         oploterror, abin.xbin, abin.meany, abin.sigymean, psym=symcat(9,thick=5), $
;           symsize=1.1, thick=6, color=fsc_color('blue',101), $
;           errcolor=fsc_color('blue',101)
;         djs_oplot, mm.xbin, mm.medy, psym=6, symsize=0.5
          legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
            '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
            /left, /top, box=0, margin=0, charsize=1.4

; local MZ          
          ohlocal = mz_closedbox(maxis,mzlocal.coeff)
          good = where(ohlocal gt yrange[0]+0.1)
          djs_oplot, maxis[good], ohlocal[good], $
            line=localline, color=im_color(localcolor), thick=6

; overplot the evolutionary model derived in FIT_MZLZEVOL
          if (iz gt 0) then begin
;            ww = where(strtrim(mzavg.calib,2) eq calib)
;            djs_oplot, maxis, mzevol_func(maxis,mzavg.mzevol_coeffs_r0zero[*,ww],$
;              z=zbins[ii].zbin,qz0=mzavg.qz0),line=2, color=im_color(calibcolor[ii],101), thick=5
;            djs_oplot, maxis, mzevol_func(maxis,mzavg.mzevol_coeffs[*,ww],$
;              z=zbins[ii].zbin,qz0=mzavg.qz0),line=2, color=im_color(calibcolor[ii],101), thick=10
             
;            ohmodel = ohlocal + (zbins[iz].zbin-mzavg.qz0)*$
;              poly(maxis-mzavg.dlogohdz_normmass,mzavg.dlogohdz_coeff)
;            keep = where((ohmodel gt yrange[0]+0.1))
;            djs_oplot, maxis[keep], ohmodel[keep], line=evolline, $
;              color=fsc_color(evolcolor,101), thick=8
          endif
          xyouts, pos[0,5], pos[1,5]-0.07, align=0.5, /norm, xtitle1
       endfor 
;      xyouts, pos[0,0]-0.075, pos[1,0], ytitle1, align=0.5, orientation=90, /norm
;      xyouts, pos[0,2]-0.075, pos[1,2], ytitle1, align=0.5, orientation=90, /norm
       im_plotconfig, /psclose, psfile=psfile
    endfor

stop    
    
; --------------------------------------------------
; compare the AGES LZ evolution to that inferred from the mock-AGES
; sample

    

return
end


;;; --------------------------------------------------
;;; Figure 18 - AGES - LZ evolution
;;    magrange1 = [-16.3,-23.7]
;;    ohrange1 = [8.35,9.37] ; [8.3,9.39]
;;
;;    localline = 0
;;    localcolor = 'black'
;;    evolline = 5
;;    evolcolor = 'firebrick'
;;    
;;    mzevol = mrdfits(mzpath+'mzevol_avg.fits.gz',1)
;;    for ii = 0, 2 do begin
;;       case ii of
;;          0: begin
;;             t04 = 1 & m91 = 0 & kk04 = 0
;;             calib = 't04'
;;             ohrange1 = [8.3,9.35]
;;          end
;;          1: begin
;;             t04 = 0 & m91 = 1 & kk04 = 0
;;             calib = 'm91'
;;             ohrange1 = [8.3,9.2]
;;          end
;;          2: begin
;;             t04 = 0 & m91 = 0 & kk04 = 1
;;             calib = 'kk04'
;;             ohrange1 = [8.55,9.35]
;;          end
;;       endcase
;;       ohtitle = mzplot_ohtitle(t04=t04,m91=m91,kk04=kk04)
;;
;;       mzlocal = mrdfits(mzpath+'mzlocal_sdss_ews_'+calib+'.fits.gz',1)
;;       lzlocal = mrdfits(mzpath+'lzlocal_sdss_ews_'+calib+'.fits.gz',1)
;;       lzevol = mrdfits(mzpath+'lzevol_'+calib+'.fits.gz',1)
;;
;;;      ainfo = mzlz_grab_info(agesohnodust,agesancillary,agesmass,$
;;;        t04=t04,m91=m91,kk04=kk04,/nolimit,/flux,zmin=0.05,zmax=0.15)
;;       ainfo = mzlz_grab_info(agesohdust,agesancillary,agesmass,$
;;         t04=t04,m91=m91,kk04=kk04,/nolimit,/errcut) ; apply an error cut for the contours
;;       
;;       psfile = pspath+'mzevol_'+calib+suffix
;;       mzplot_sixpanel, ainfo.z, ainfo.mass, ainfo.oh, ainfo.weight, $
;;         oh_err=ainfo.oh_err, psfile=psfile, xtitle=mzplot_masstitle(), $
;;         ytitle=ohtitle, /ages, xrange=massrange1, yrange=ohrange1, npix=10, $
;;         mzlocal=mzlocal, mzevol=mzevol, localline=localline, $
;;         localcolor=localcolor, evolline=evolline, evolcolor=evolcolor, $
;;         postscript=keyword_set(ps)
;;       
;;       psfile = pspath+'lzevol_'+calib+suffix
;;       mzplot_sixpanel, ainfo.z, ainfo.mb_ab, ainfo.oh, ainfo.weight, $
;;         oh_err=ainfo.oh_err, psfile=psfile, xtitle=mzplot_mbtitle(), $
;;         ytitle=ohtitle, /ages, xrange=magrange1, yrange=ohrange1, npix=10, $
;;         lzlocal=lzlocal, lzevol=lzevol, localline=localline, $
;;         localcolor=localcolor, evolline=evolline, evolcolor=evolcolor, $
;;         postscript=keyword_set(ps)
;;    endfor
;;       
