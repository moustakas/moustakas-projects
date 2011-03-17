pro write_deep2_oiilf_sample, kcorr, ispec, write=write
; jm07sep27nyu - written based on WRITE_AGES_MZ_SAMPLE
; jm08may02nyu - updated
    
    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33

    oiilfpath = deep2_path(/projects)+'oiilf/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = oiilfpath+'write_deep2_oiilf_sample.log'
       splog, filename=splogfile
       splog, 'Log file '+splogfile+' opened '+systime()
    endif
    splog, 'IDL version: ' + string(!version,format='(99(A," "))')
    spawn, 'uname -a', uname
    splog, 'UNAME: '+uname[0]
       
    if keyword_set(write) then begin
       postscript = 1L
       postthick1 = 4.0
       postthick2 = 8.0
    endif else begin
       postthick1 = 1.8
       postthick2 = 2.0
       im_window, 0, xratio=0.45, /square
    endelse

    charsize_6 = 1.6

    ewoiicut = 2.0 ; 10.0
    snrcut1 = 10.0 ; 5.0 ; for sample selection
    snrcut2 = 10.0 ; for plot testing, below

    ewaxis = findgen((alog10(1E3)-alog10(1E-3))/0.05)*0.05+alog10(1E-3)
    ratiorange = [0.3,5.0]
    snrrange = [0.02,300.0]
    ewrange = [0.03,1000.0]

; ---------------------------------------------------------------------------    
; read the data
; ---------------------------------------------------------------------------    

    splog, 'Reading the data...'
    if (n_elements(kcorr) eq 0L) then kcorr = read_deep2(/kcorr)
    if (n_elements(ispec) eq 0L) then begin
       ispec1 = read_deep2(/specdata)
       ispecmore = {oii_lum: -999.0, oii_lum_err: -999.0}
       ispecmore = replicate(ispecmore,n_elements(ispec1))
       ispec = struct_addtags(ispec1,ispecmore)
    endif
    ngalaxy = n_elements(kcorr)

; ---------------------------------------------------------------------------
; read the VMAX values
; ---------------------------------------------------------------------------

;   fattempt = mrdfits(deep2_path(/analysis)+'deep2_completeness.fits',0,/silent)
;   fgot = mrdfits(deep2_path(/analysis)+'deep2_completeness.fits',1,/silent)
;   num = 15L
;   rmilim = [0.0,2.0]
;   rlim = [19.0,24.1]

; ---------------------------------------------------------------------------
; parent sample cuts
; ---------------------------------------------------------------------------

    survey_area = 3.0*!dtor^2.0 ; 3.05D-4 sr/deg^2 [sr]
    survey_filter = 'deep_R.par'
    sample_zmin = 0.7
    sample_zmax = 1.5
    deep2_mlimit = [19.0,24.1]
    
    parent = where((kcorr.z gt sample_zmin) and (kcorr.z lt sample_zmax),nparent)
    splog, 'Parent sample: '+string(nparent,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      strtrim(string(100.0*nparent/float(ngalaxy),format='(F12.1)'),2)+'%).'

    parent_kcorr = kcorr[parent]
    parent_ispec = ispec[parent]

; statistics    
    
    zstats = im_stats(parent_kcorr.z)
    splog, '   Redshift      : ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; ---------------------------------------------------------------------------
; define the [OII] sample
; ---------------------------------------------------------------------------
    
    oiilf = where($
      (parent_ispec.oii_3727_1[0]/parent_ispec.oii_3727_1[1] gt snrcut1) and $
      (parent_ispec.oii_3727_2[0]/parent_ispec.oii_3727_2[1] gt snrcut1) and $
      (parent_ispec.oii_3727_1_ew[0] gt ewoiicut) and (parent_ispec.oii_3727_2_ew[0] gt ewoiicut),noiilf)
    splog, '[OII] sample: '+string(noiilf,format='(I0)')+'/'+string(nparent,format='(I0)')+' ('+$
      strtrim(string(100.0*noiilf/nparent,format='(F12.1)'),2)+'%).'

    oiilf_kcorr = parent_kcorr[oiilf]
    oiilf_ispec = parent_ispec[oiilf]

; statistics    
    
    zstats = im_stats(oiilf_kcorr.z)
    splog, '   Redshift      : ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

;; ---------------------------------------------------------------------------
;; compute the statistical weight for each galaxy
;; ---------------------------------------------------------------------------
;
;    rmi = deep2dust[parent].magr-deep2dust[parent].magi
;    rmag = deep2dust[parent].magr
;
;    irmi = long((rmi-rmilim[0])/(rmilim[1]-rmilim[0])*float(num))
;    ir = long((rmag-rlim[0])/(rlim[1]-rlim[0])*float(num))
;    crap = where((irmi ge num) or (irmi lt 0) or (ir ge num) or (ir lt 0),comp=good,ncomp=ngood)
;    if (ngood ne 0L) then moretags[good].weight = fgot[ir[good],irmi[good]]*fattempt[ir[good],irmi[good]]
    
; ---------------------------------------------------------------------------
; compute the line-luminosity 
; ---------------------------------------------------------------------------

    splog, 'Computing emission-line luminosities'

    dlum = dluminosity(oiilf_kcorr.z,/cm)
    oii_lum_factor = oiilf_kcorr.cflux_3727*(4.0*!dpi*dlum^2.0)/lsun

    oiilf_ispec.oii_lum = alog10(oiilf_ispec.oii_3727_ew[0]*oii_lum_factor)
    oiilf_ispec.oii_lum_err = oiilf_ispec.oii_3727_ew[1]/oiilf_ispec.oii_3727_ew[0]/alog(10.0)

; ---------------------------------------------------------------------------
; sample selection plots
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin
       psname = oiilfpath+'deep2_oiilf_sample_selection.ps'
       dfpsplot, psname, /color, /square
    endif

; some useful variables    
    
    ewoii_1 = parent_ispec.oii_3727_1_ew[0] & ewoii_1_err = parent_ispec.oii_3727_1_ew[1]
    snroii_1 = parent_ispec.oii_3727_1[0]/parent_ispec.oii_3727_1[1] ; note FLUX, not EW!
    oiigood_1 = where((ewoii_1 gt 0.0) and (snroii_1 gt 0.0))

    ewoii_2 = parent_ispec.oii_3727_2_ew[0] & ewoii_2_err = parent_ispec.oii_3727_2_ew[1]
    snroii_2 = parent_ispec.oii_3727_2[0]/parent_ispec.oii_3727_2[1] ; note FLUX, not EW!
    oiigood_2 = where((ewoii_2 gt 0.0) and (snroii_2 gt 0.0))

; ------------------------------------------------------------
; EW([O II] 3726) vs EW([O II] 3729), coded by S/N

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.6, height=3.3*[1,1], $
      xmargin=[1.6,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.1, $
      position=pos, /normal

    xtitle = 'EW([O II] \lambda3729)'
    ytitle = 'EW([O II] \lambda3726)'
    xrange = ewrange
    yrange = xrange

    hi = where((snroii_1 gt snrcut2) and (snroii_2 gt snrcut2),comp=lo)

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog, $
      xtickname=replicate(' ',10)
    djs_oplot, 10^ewaxis, poly(10^ewaxis,[0.0,0.67]), line=0, thick=postthick1 ; ~1 cm^-3
    djs_oplot, 10^ewaxis, poly(10^ewaxis,[0.0,2.32]), line=5, thick=postthick1 ; ~5000 cm^-3
    djs_oplot, ewoii_2[hi], ewoii_1[hi], ps=3, color='green'
    legend, 'S/N > '+string(snrcut2,format='(I0)'), /left, /top, box=0, charsize=1.5, $
      charthick=postthick1, psym=3, color=djs_icolor('green')

    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1], /xlog, /ylog
    djs_oplot, 10^ewaxis, poly(10^ewaxis,[0.0,0.67]), line=0, thick=postthick1 ; ~1 cm^-3
    djs_oplot, 10^ewaxis, poly(10^ewaxis,[0.0,2.32]), line=5, thick=postthick1 ; ~5000 cm^-3
    djs_oplot, ewoii_2[lo], ewoii_1[lo], ps=3, color='cyan'
    legend, 'S/N < '+string(snrcut2,format='(I0)'), /left, /top, box=0, charsize=1.5, $
      charthick=postthick1, psym=3, color=djs_icolor('cyan')
    
    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; EW([O II]) vs S/N EW([OII])

    trial_snrcut = [3.0,5.0,10.0,20.0]
    
; [O II] 3726

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.6, height=3.3*[1,1], $
      xmargin=[1.6,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.1, $
      position=pos, /normal

    xtitle = 'EW([O II]) (\AA)'
    ytitle = 'S/N EW([O II] \lambda3726)'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog, $
      xtickname=replicate(' ',10)

    oiicoeff_1 = robust_linefit(alog10(ewoii_1[oiigood_1]),$
      alog10(snroii_1[oiigood_1]),/bisect,oiiyfit_1,oiisig_1)
    oiiaxis_1 = poly(ewaxis,oiicoeff_1)
    djs_oplot, ewoii_1, snroii_1, psym=3, color='blue'
    djs_oplot, 10.0^ewaxis, 10.0^oiiaxis_1, line=0, thick=postthick2

    plot_ewoiicut_1 = interpol(10.0^ewaxis,10.0^oiiaxis_1,trial_snrcut)
    splog, 'EW([O II] 3726) at S/N = '+strjoin(string(trial_snrcut,format='(I0)'),', ')+': ', $
      plot_ewoiicut_1

    for it = 0L, n_elements(trial_snrcut)-1L do begin
       djs_oplot, plot_ewoiicut_1[it]*[1,1], [10^!y.crange[0],trial_snrcut[it]], line=2, thick=postthick2
       djs_oplot, [10^!x.crange[0],plot_ewoiicut_1[it]], trial_snrcut[it]*[1,1], line=2, thick=postthick2
    endfor
       
; [O II] 3729

    ytitle = 'S/N EW([O II] \lambda3729)'

    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1], /xlog, /ylog

    oiicoeff_2 = robust_linefit(alog10(ewoii_2[oiigood_2]),$
      alog10(snroii_2[oiigood_2]),/bisect,oiiyfit_2,oiisig_2)
    oiiaxis_2 = poly(ewaxis,oiicoeff_2)
    djs_oplot, ewoii_2, snroii_2, psym=3, color='red'
    djs_oplot, 10.0^ewaxis, 10.0^oiiaxis_2, line=0, thick=postthick2

    plot_ewoiicut_2 = interpol(10.0^ewaxis,10.0^oiiaxis_2,trial_snrcut)
    splog, 'EW([O II] 3729) at S/N = '+strjoin(string(trial_snrcut,format='(I0)'),', ')+': ', $
      plot_ewoiicut_2

    for it = 0L, n_elements(trial_snrcut)-1L do begin
       djs_oplot, plot_ewoiicut_2[it]*[1,1], [10^!y.crange[0],trial_snrcut[it]], line=2, thick=postthick2
       djs_oplot, [10^!x.crange[0],plot_ewoiicut_2[it]], trial_snrcut[it]*[1,1], line=2, thick=postthick2
    endfor
       
;   legend, textoidl('EW([O II]) > '+string(ewoiicut,format='(G0.0)')+' \AA'), /left, /top, $
;     box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; redshift vs EW([O II]) and [O II] doublet ratio

    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.6, height=3.3*[1,1], $
      xmargin=[1.6,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.1, $
      position=pos, /normal

    xtitle = 'Redshift'
    ytitle1 = 'EW([O II]) [\AA]'
    ytitle2 = '[O II] \lambda3726 / [O II] \lambda3729'
    xrange = [sample_zmin,sample_zmax]
    yrange1 = ewrange
    yrange2 = ratiorange

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange1, $
      xtitle='', ytitle=ytitle1, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      xtickname=replicate(' ',10), /ylog
    djs_oplot, parent_kcorr.z, parent_ispec.oii_3727_ew[0], psym=3, color='grey'
    djs_oplot, oiilf_kcorr.z, oiilf_ispec.oii_3727_ew[0], psym=3, color='dark green'

    djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, xrange=xrange, yrange=yrange2, $
      xtitle=xtitle, ytitle=ytitle2, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,1], /ylog
    djs_oplot, oiilf_kcorr.z, 1.05*oiilf_ispec.oii_3727_1_ew[0]/oiilf_ispec.oii_3727_2_ew[0], $
      psym=3, color='dark green'
    djs_oplot, !x.crange, 0.67*[1,1], line=0, thick=2 ; ~1 cm^-3
    djs_oplot, !x.crange, 2.32*[1,1], line=5, thick=2 ; ~5000 cm^-3
    
    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; close the PS file    
    
    if keyword_set(write) then begin
       dfpsclose
       spawn, 'gzip -f '+psname
    endif

; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin

; K-corrections       
       
       parentfile_kcorr = oiilfpath+'deep2_oiilf_parent_kcorr.fits'
       splog, 'Writing '+parentfile_kcorr
       mwrfits, parent_kcorr, parentfile_kcorr, /create

       oiilffile_kcorr = oiilfpath+'deep2_oiilf_kcorr.fits'
       splog, 'Writing '+oiilffile_kcorr
       mwrfits, oiilf_kcorr, oiilffile_kcorr, /create

; emission-lilf data

       oiilffile = oiilfpath+'deep2_oiilf_ispec.fits'
       splog, 'Writing '+oiilffile
       mwrfits, struct_addtags(oiilf_temden,oiilf_ispec), oiilffile, /create

       splog, 'Total time = '+strtrim(string((systime(1)-$
         stime0)/60.0,format='(F12.1)'),2)+' minutes.'
       splog, /close

    endif

stop
    
return
end
