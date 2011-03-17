pro write_sdss_oiilf_sample, kcorr, ispec, postlss=postlss, dr6=dr6, write=write
; jm07sep27nyu - written based on WRITE_AGES_MZ_SAMPLE
; jm08may02nyu - updated
    
    q0 = 1.50 & q1 = 0.0 & qz0 = 0.1

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33

    oiilfpath = deep2_path(/projects)+'oiilf/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = oiilfpath+'write_sdss_oiilf_sample.log'
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

    ewoiicut = 4.5 ; 10.0
    snrcut1 = 5.0

    ewaxis = findgen((alog10(1E3)-alog10(1E-3))/0.05)*0.05+alog10(1E-3)
    ratiorange = [0.3,5.0]
    snrrange = [0.02,300.0]
    ewrange = [0.03,1000.0]

; ---------------------------------------------------------------------------    
; read the data
; ---------------------------------------------------------------------------    

    splog, 'Reading the data...'
    if (n_elements(kcorr) eq 0L) then begin
       kcorr1 = read_sdss_main(/kcorr,dr6=dr6)
       kcorr = struct_addtags(kcorr1,replicate({vol: -999.0, $
         zmin_noevol: -999.0, zmax_noevol: -999.0, vmax_noevol: -999.0, $
         zmin_levol: -999.0, zmax_levol: -999.0, vmax_levol: -999.0, $
         object_position: 0L},n_elements(kcorr1)))
    endif
    if (n_elements(ispec) eq 0L) then begin
       ispec1 = read_sdss_main(/ispec,dr6=dr6)
       ispecmore = {oii_lum: -999.0, oii_lum_err: -999.0}
       ispecmore = replicate(ispecmore,n_elements(ispec1))
       ispec = struct_addtags(ispec1,ispecmore)
    endif
    ngalaxy = n_elements(kcorr)

; we need OBJECT_POSITION from POSTLSS, so copy that over right now
    
    if (n_elements(postlss) eq 0L) then postlss = read_sdss_main(/postlss,dr6=dr6)
    kcorr.object_position = postlss.object_position

; ---------------------------------------------------------------------------
; parent sample cuts
; ---------------------------------------------------------------------------

    if keyword_set(dr6) then dr_area = 6750.0*!dtor^2.0 else $
      dr_area = 4783.0*!dtor^2.0 ; 3.05D-4 sr/deg^2 [sr]
    vname = 'default.nolines' ; this should match WRITE_SDSS_MAIN
    survey_filter = 'sdss_r0.par'
    sample_zmin = 0.04  ; ensure [O II] is in range
    sample_zmax = 0.30  ; very few objects above this redshift
    sdss_mlimit = [14.5,17.6]
    sdss_band_shift = 0.1
    sdss_r_aboffset = 0.01
    
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
; compute zmin, zmax, and vmax, but only when writing out
; ---------------------------------------------------------------------------    

    if keyword_set(write) then begin

       parent_kcorr.vol  = (dr_area/3.0)*(lf_comvol(parent_kcorr.z)-$
         (lf_comvol(sample_zmin))[0])/h100^3.0 ; h=1-->h=0.7 [Mpc^3]

       splog, 'Computing Vmax - no luminosity evolution'
       lf_calc_vmax, -2.5*alog10(parent_kcorr.abmaggies[2]), parent_kcorr.ugriz_absmag[2]-5.0*alog10(h100), $ ; h=0.7-->h=1
         parent_kcorr.coeffs, survey_filter, dr_area, sdss_mlimit[0], sdss_mlimit[1], $
         sample_zmin, sample_zmax, vname=vname, band_shift=sdss_band_shift, magoffset=sdss_r_aboffset, $
         actual_z=parent_kcorr.z, vmax=vmax, zmin=zmin, zmax=zmax

       parent_kcorr.zmin_noevol = zmin
       parent_kcorr.zmax_noevol = zmax
       parent_kcorr.vmax_noevol = vmax/h100^3.0 ; h=1-->h=0.7 [Mpc^3]

       splog, 'Computing Vmax - with luminosity evolution'
       lf_calc_vmax, -2.5*alog10(parent_kcorr.abmaggies[2]), parent_kcorr.ugriz_absmag[2]-5.0*alog10(h100), $ ; h=0.7-->h=1
         parent_kcorr.coeffs, survey_filter, dr_area, sdss_mlimit[0], sdss_mlimit[1], $
         sample_zmin, sample_zmax, vname=vname, band_shift=sdss_band_shift, magoffset=sdss_r_aboffset, $
         actual_z=parent_kcorr.z, vmax=vmax, zmin=zmin, zmax=zmax, $
         q0=q0, q1=q1, qz0=qz0

       parent_kcorr.zmin_levol = zmin
       parent_kcorr.zmax_levol = zmax
       parent_kcorr.vmax_levol = vmax/h100^3.0 ; h=1-->h=0.7 [Mpc^3]

    endif

; ---------------------------------------------------------------------------
; define the [OII] sample
; ---------------------------------------------------------------------------
    
    oiilf = where($
      (parent_ispec.oii_3727[0]/parent_ispec.oii_3727[1] gt snrcut1) and $
      (parent_ispec.oii_3727_ew[0] gt ewoiicut),noiilf)
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
;    rmi = sdssdust[parent].magr-sdssdust[parent].magi
;    rmag = sdssdust[parent].magr
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
       psname = oiilfpath+'sdss_oiilf_sample_selection.ps'
       dfpsplot, psname, /color, /square
    endif

; some useful variables    
    
    ewoii = parent_ispec.oii_3727_ew[0] & ewoii_err = parent_ispec.oii_3727_ew[1]
    snroii = parent_ispec.oii_3727[0]/parent_ispec.oii_3727[1] ; note FLUX, not EW!
    oiigood = where((ewoii gt 0.0) and (snroii gt 0.0))

; ------------------------------------------------------------
; EW([O II]) vs S/N EW([OII])

    trial_snrcut = [3.0,5.0,10.0,20.0]
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.6,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.1, $
      position=pos, /normal

    xtitle = 'EW([O II]) (\AA)'
    ytitle = 'S/N EW([O II])'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog

    oiicoeff = robust_linefit(alog10(ewoii[oiigood]),$
      alog10(snroii[oiigood]),/bisect,oiiyfit,oiisig)
    oiiaxis = poly(ewaxis,oiicoeff)
    djs_oplot, ewoii, snroii, psym=3, color='red'
    djs_oplot, 10.0^ewaxis, 10.0^oiiaxis, line=0, thick=postthick2

    plot_ewoiicut = interpol(10.0^ewaxis,10.0^oiiaxis,trial_snrcut)
    splog, 'EW([O II]) at S/N = '+strjoin(string(trial_snrcut,format='(I0)'),', ')+': ', $
      plot_ewoiicut

    for it = 0L, n_elements(trial_snrcut)-1L do begin
       djs_oplot, plot_ewoiicut[it]*[1,1], [10^!y.crange[0],trial_snrcut[it]], line=2, thick=postthick2
       djs_oplot, [10^!x.crange[0],plot_ewoiicut[it]], trial_snrcut[it]*[1,1], line=2, thick=postthick2
    endfor
       
;   legend, textoidl('EW([O II]) > '+string(ewoiicut,format='(G0.0)')+' \AA'), /left, /top, $
;     box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

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
       
       parentfile_kcorr = oiilfpath+'sdss_oiilf_parent_kcorr.fits'
       splog, 'Writing '+parentfile_kcorr
       mwrfits, parent_kcorr, parentfile_kcorr, /create

       oiilffile_kcorr = oiilfpath+'sdss_oiilf_kcorr.fits'
       splog, 'Writing '+oiilffile_kcorr
       mwrfits, oiilf_kcorr, oiilffile_kcorr, /create

; emission-line data

       oiilffile = oiilfpath+'sdss_oiilf_ispec.fits'
       splog, 'Writing '+oiilffile
       mwrfits, struct_addtags(oiilf_temden,oiilf_ispec), oiilffile, /create

       splog, 'Total time = '+strtrim(string((systime(1)-$
         stime0)/60.0,format='(F12.1)'),2)+' minutes.'
       splog, /close

    endif

stop
    
return
end
