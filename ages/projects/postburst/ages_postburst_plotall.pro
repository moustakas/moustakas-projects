pro ages_postburst_plotall, spec
; jm09feb10nyu - plot the final sample in one plot

    path = ages_path(/projects)+'postburst/brown/'
    sum = rsex(path+'summary_09feb23.txt')
;   sum = rsex(path+'summary.txt')

    aa1 = read_ages(/ancillary)
    ii1 = read_ages(/ispec)

    m1 = lonarr(n_elements(sum))
    for kk = 0L, n_elements(sum)-1L do m1[kk] = where(aa1.ages_id eq sum[kk].ages_id)
    aa = aa1[m1]
    ii = ii1[m1]
    absmag = sum.absmag
    
;   spherematch, aa1.ra, aa1.dec, 15.0*im_hms2dec(sum.ra), $
;     im_hms2dec(sum.dec), 1.0/3600.0, m1, m2
;
;   aa = aa1[m1]
;   ii = ii1[m1]
;   absmag = sum[m2].absmag
;   srt = sort(absmag)
;   aa = aa[srt]
;   ii = ii[srt]
;   sum = sum[m2[srt]]
;   absmag = absmag[srt]

    if (n_elements(spec) eq 0L) then $
      spec = read_ages_specfit(aa.galaxy)
    sz = size(spec,/dim)
    
; 6x4 panels       

    ncols = 6L & nrows = 4L
    xmargin = [0.9,0.2] & ymargin = [0.2,1.0]
    xspace = 0.0 & yspace = 0.0
    width = 1.65
    height = 1.825
    xpage = 8.5 & ypage = 11.0
    thislandscape = 1

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage

    psfile = path+'postburst_plotall_09feb23.ps'
    splog, 'Writing '+psfile
    im_plotfaves, /post, charthick=4.0

    arm_plotconfig, landscape=thislandscape, nx=ncols, ny=nrows, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L;, /show
;   cleanplot, /silent

    xtitle = textoidl('Observed Wavelength (\AA)')
    ytitle = 'Relative Flux (arbitrary scale)'
    charsize1 = 1.6

    xrange = [3800,8600]
    yrange = [0.01,2.5]

    for ii = 0L, sz[2]-1L do begin

; now make the plot       
       
       if (ii mod ncols) eq 0L then begin
          delvarx, ytickname
          ytitle1 = ytitle
       endif else begin
          ytickname = replicate(' ',10)
          ytitle1 = ''
       endelse

       if (ii ge 18) then begin
          delvarx, xtickname
       endif else begin
          xtickname = replicate(' ',10)
       endelse

;      npix = (size(spec,/dim))[0]
;      wave = rebin(spec[*,0,ii],npix/3)
;      flux = rebin(spec[*,1,ii],npix/3)
;      modelflux = rebin(spec[*,2,ii]+spec[*,4,ii],npix/3)
;      good = where(wave gt 0.0)
;      norm = median(flux[good])
;      modelflux = modelflux[good]/norm
;      flux = flux[good]/norm

       good = where(spec[*,0,ii] gt 0.0,npix)
       wave = spec[good,0,ii]*(1.0+aa[ii].z)
       flux = spec[good,1,ii]-spec[good,4,ii]
       modelflux = spec[good,2,ii]
;      flux = spec[good,1,ii]
;      modelflux = spec[good,2,ii]+spec[good,4,ii]
       norm = median(flux)
       modelflux = modelflux/norm
       flux = flux/norm
;      yrange = im_max(flux,sigrej=3.0)*[-0.1,1.5]

;      bfactor = 4
;      if odd(npix) then begin
;         wave = rebin(wave[0:npix-2],npix/bfactor)
;         flux = rebin(flux[0:npix-2],npix/bfactor)
;         modelflux = rebin(modelflux[0:npix-2],npix/bfactor)
;      endif else begin
;         wave = rebin(wave,npix/bfactor)
;         flux = rebin(flux,npix/bfactor)
;         modelflux = rebin(modelflux,npix/bfactor)
;      endelse

       bfactor = 4
       chop = npix mod 1024
       wave = rebin(wave[0:npix-chop-1],(npix-chop)/bfactor)
       flux = rebin(flux[0:npix-chop-1],(npix-chop)/bfactor)
       modelflux = rebin(modelflux[0:npix-chop-1],(npix-chop)/bfactor)

       djs_plot, [0], [0], /nodata, charsize=charsize1, xtitle='', ytitle='', $
         xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,ii], $
         ytickname=ytickname, xtickname=xtickname, noerase=(ii gt 0L), $
         xtickinterval=2500, ytickinterval=1.0
;      legend, strtrim(string(absmag[ii],format='(F12.1)'),2), /right, /top, $
;        box=0, charsize=1.0, /clear, margin=0
       case ii of
          0: begin
             up = 0.4
             dn = -0.3
          end
          1: begin
             up = 0.3
             dn = -0.4
          end
          2: begin
             up = 0.4
             dn = -0.3
          end
          5: begin
             up = 0.5
             dn = -0.2
          end
          23: begin
             up = 0.5
             dn = -0.2
          end
;         21: begin
;            up = 0.3
;            dn = -0.2
;         end
;         22: begin
;            up = 0.1
;            dn = -0.5
;         end
;         23: begin
;            up = 0.3
;            dn = -0.2
;         end
          else: begin
             up = 0.6
             dn = -0.1
          end
       endcase
       djs_oplot, wave, smooth(flux,5)+up, ps=10, thick=1.5
       djs_oplot, wave, smooth(modelflux,5)+dn, ps=10, thick=1.5, color='red'

    endfor

; x-title    
;   xyouts, pos[0,19], pos[1,19]-0.06, xtitle, align=0.5, charsize=1.4, /normal
    xyouts, pos[0,21], pos[1,21]-0.08, xtitle, align=0.5, charsize=charsize1, /normal
;   xyouts, pos[0,23], pos[1,23]-0.06, xtitle, align=0.5, charsize=1.4, /normal
; y-title    
;   xyouts, pos[0,0]-0.05, pos[1,0], ytitle, align=0.5, orientation=90, charsize=1.4, /normal
    xyouts, pos[0,6]-0.04, pos[1,6], ytitle, align=0.5, orientation=90, charsize=charsize1, /normal

    dfpsclose
    im_plotfaves

stop    
    
return
end
    
