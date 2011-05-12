pro bcgs_plots
; jm10jul22ucsd - basic plots for the BCGS project with Anthony
; jm11apr07ucsd - major updates

    bcgspath = ages_path(/projects)+'bcgs/'
    isedpath = bcgspath+'isedfit/'
    qapath = bcgspath+'qaplots/'

; read all the isedfit output files
    bc03_solar = read_bcgs_isedfit(supergrid=1)
    basti_ss_solar = read_bcgs_isedfit(supergrid=3)
    basti_ae_solar = read_bcgs_isedfit(supergrid=5)
    
    bc03_supersolar = read_bcgs_isedfit(supergrid=2)
    basti_ss_supersolar = read_bcgs_isedfit(supergrid=4)
    basti_ae_supersolar = read_bcgs_isedfit(supergrid=6)
    
    sdss_bc03_solar = read_bcgs_isedfit(supergrid=1,/sdss)
    sdss_basti_ss_solar = read_bcgs_isedfit(supergrid=3,/sdss)
    sdss_basti_ae_solar = read_bcgs_isedfit(supergrid=5,/sdss)
    
    sdss_bc03_supersolar = read_bcgs_isedfit(supergrid=2,/sdss)
    sdss_basti_ss_supersolar = read_bcgs_isedfit(supergrid=4,/sdss)
    sdss_basti_ae_supersolar = read_bcgs_isedfit(supergrid=6,/sdss)

; plotting preferences    
    bc03_color = 'dark green'      & bc03_psym = 6
    basti_ss_color = 'dodger blue' & basti_ss_psym = 9
    basti_ae_color = 'firebrick'   & basti_ae_psym = 5

    sdss_bc03_color = 'dark green'      & sdss_bc03_psym = 15
    sdss_basti_ss_color = 'dodger blue' & sdss_basti_ss_psym = 16
    sdss_basti_ae_color = 'firebrick'   & sdss_basti_ae_psym = 17
    
; --------------------------------------------------
; age vs redshift
    zrange = [-0.05,2]
    agerange = [0.0,14.0]
    ztitle = 'Redshift'
    agetitle = 'Age (Gyr)'
    zformmax = 5.0
    symsize = 1.2

    zaxis = range(0.0,2.0,100)

    showerr = 0
    
    psfile = qapath+'z_vs_age.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], $
      width=6.8

; -------------------------=
; solar
    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=agerange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=agetitle

; bc03
    if keyword_set(showerr) then begin
       oploterror, sdss_bc03_solar.zobj, sdss_bc03_solar.age_50, sdss_bc03_solar.age_err, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100), $
         errcolor=fsc_color(sdss_bc03_color,100), symsize=symsize1, errthick=4
       oploterror, bc03_solar.zobj, bc03_solar.age_50, bc03_solar.age_err, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100), $
         errcolor=fsc_color(bc03_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_bc03_solar.zobj, sdss_bc03_solar.age_50, symsize=symsize1, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100)
       oplot, bc03_solar.zobj, bc03_solar.age_50, symsize=symsize1, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100)
    endelse

; basti-ss
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ss_solar.zobj, sdss_basti_ss_solar.age_50, sdss_basti_ss_solar.age_err, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100), $
         errcolor=fsc_color(sdss_basti_ss_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ss_solar.zobj, basti_ss_solar.age_50, basti_ss_solar.age_err, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100), $
         errcolor=fsc_color(basti_ss_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ss_solar.zobj, sdss_basti_ss_solar.age_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100)
       oplot, basti_ss_solar.zobj, basti_ss_solar.age_50, symsize=symsize1, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100)
    endelse

; basti-ae
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ae_solar.zobj, sdss_basti_ae_solar.age_50, sdss_basti_ae_solar.age_err, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100), $
         errcolor=fsc_color(sdss_basti_ae_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ae_solar.zobj, basti_ae_solar.age_50, basti_ae_solar.age_err, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100), $
         errcolor=fsc_color(basti_ae_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ae_solar.zobj, sdss_basti_ae_solar.age_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100)
       oplot, basti_ae_solar.zobj, basti_ae_solar.age_50, symsize=symsize1, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100)
    endelse

    im_legend, ['BC03','BaSTI-ss','BaSTI-ae'], /right, /top, box=0, $
      color=[bc03_color,basti_ss_color,basti_ae_color], $
      psym=[bc03_psym,basti_ss_psym,basti_ae_psym], $
      symthick=4, charsize=1.5
    xyouts, 0.6, 12.0, 'Z = 0.02', align=0.5, /data, charsize=1.5

    djs_oplot, zaxis, getage(zaxis), line=0, thick=6
    djs_oplot, zaxis, getage(zaxis)-getage(zformmax), line=5, thick=6

; -------------------------=
; super-solar
    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=agerange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=agetitle

; bc03
    if keyword_set(showerr) then begin
       oploterror, sdss_bc03_supersolar.zobj, sdss_bc03_supersolar.age_50, sdss_bc03_supersolar.age_err, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100), $
         errcolor=fsc_color(sdss_bc03_color,100), symsize=symsize1, errthick=4
       oploterror, bc03_supersolar.zobj, bc03_supersolar.age_50, bc03_supersolar.age_err, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100), $
         errcolor=fsc_color(bc03_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_bc03_supersolar.zobj, sdss_bc03_supersolar.age_50, symsize=symsize1, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100)
       oplot, bc03_supersolar.zobj, bc03_supersolar.age_50, symsize=symsize1, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100)
    endelse

; basti-ss
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ss_supersolar.zobj, sdss_basti_ss_supersolar.age_50, sdss_basti_ss_supersolar.age_err, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100), $
         errcolor=fsc_color(sdss_basti_ss_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ss_supersolar.zobj, basti_ss_supersolar.age_50, basti_ss_supersolar.age_err, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100), $
         errcolor=fsc_color(basti_ss_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ss_supersolar.zobj, sdss_basti_ss_supersolar.age_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100)
       oplot, basti_ss_supersolar.zobj, basti_ss_supersolar.age_50, symsize=symsize1, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100)
    endelse

; basti-ae
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ae_supersolar.zobj, sdss_basti_ae_supersolar.age_50, sdss_basti_ae_supersolar.age_err, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100), $
         errcolor=fsc_color(sdss_basti_ae_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ae_supersolar.zobj, basti_ae_supersolar.age_50, basti_ae_supersolar.age_err, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100), $
         errcolor=fsc_color(basti_ae_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ae_supersolar.zobj, sdss_basti_ae_supersolar.age_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100)
       oplot, basti_ae_supersolar.zobj, basti_ae_supersolar.age_50, symsize=symsize1, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100)
    endelse

    im_legend, ['BC03','BaSTI-ss','BaSTI-ae'], /right, /top, box=0, $
      color=[bc03_color,basti_ss_color,basti_ae_color], $
      psym=[bc03_psym,basti_ss_psym,basti_ae_psym], $
      symthick=4, charsize=1.5
    xyouts, 0.6, 12.0, 'Z = 0.05', align=0.5, /data, charsize=1.5

    djs_oplot, zaxis, getage(zaxis), line=0, thick=6
    djs_oplot, zaxis, getage(zaxis)-getage(zformmax), line=5, thick=6

    im_plotconfig, psfile=psfile, /psclose, /gzip


; --------------------------------------------------
; mass vs redshift
    zrange = [-0.05,2]
    massrange = [10.0,12.6]
    ztitle = 'Redshift'
    masstitle = 'log (M/M_{'+sunsymbol()+'})'
    symsize = 1.2

    zaxis = range(0.0,2.0,100)
    showerr = 0
    
    psfile = qapath+'z_vs_mass.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], $
      width=6.8

; -------------------------=
; solar
    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=masstitle

; bc03
    if keyword_set(showerr) then begin
       oploterror, sdss_bc03_solar.zobj, sdss_bc03_solar.mass_50, sdss_bc03_solar.mass_err, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100), $
         errcolor=fsc_color(sdss_bc03_color,100), symsize=symsize1, errthick=4
       oploterror, bc03_solar.zobj, bc03_solar.mass_50, bc03_solar.mass_err, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100), $
         errcolor=fsc_color(bc03_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_bc03_solar.zobj, sdss_bc03_solar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100)
       oplot, bc03_solar.zobj, bc03_solar.mass_50, symsize=symsize1, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100)
    endelse

; basti-ss
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ss_solar.zobj, sdss_basti_ss_solar.mass_50, sdss_basti_ss_solar.mass_err, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100), $
         errcolor=fsc_color(sdss_basti_ss_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ss_solar.zobj, basti_ss_solar.mass_50, basti_ss_solar.mass_err, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100), $
         errcolor=fsc_color(basti_ss_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ss_solar.zobj, sdss_basti_ss_solar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100)
       oplot, basti_ss_solar.zobj, basti_ss_solar.mass_50, symsize=symsize1, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100)
    endelse

; basti-ae
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ae_solar.zobj, sdss_basti_ae_solar.mass_50, sdss_basti_ae_solar.mass_err, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100), $
         errcolor=fsc_color(sdss_basti_ae_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ae_solar.zobj, basti_ae_solar.mass_50, basti_ae_solar.mass_err, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100), $
         errcolor=fsc_color(basti_ae_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ae_solar.zobj, sdss_basti_ae_solar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100)
       oplot, basti_ae_solar.zobj, basti_ae_solar.mass_50, symsize=symsize1, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100)
    endelse

    im_legend, ['BC03','BaSTI-ss','BaSTI-ae'], /right, /top, box=0, $
      color=[bc03_color,basti_ss_color,basti_ae_color], $
      psym=[bc03_psym,basti_ss_psym,basti_ae_psym], $
      symthick=4, charsize=1.5
    xyouts, 0.8, 12.2, 'Z = 0.02', align=0.5, /data, charsize=1.5

; -------------------------=
; super-solar
    djs_plot, [0], [0], /nodata, position=pos, $
      xrange=zrange, yrange=massrange, xsty=1, ysty=1, $
      xtitle=ztitle, ytitle=masstitle

; bc03
    if keyword_set(showerr) then begin
       oploterror, sdss_bc03_supersolar.zobj, sdss_bc03_supersolar.mass_50, sdss_bc03_supersolar.mass_err, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100), $
         errcolor=fsc_color(sdss_bc03_color,100), symsize=symsize1, errthick=4
       oploterror, bc03_supersolar.zobj, bc03_supersolar.mass_50, bc03_supersolar.mass_err, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100), $
         errcolor=fsc_color(bc03_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_bc03_supersolar.zobj, sdss_bc03_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_bc03_psym,thick=6), color=fsc_color(sdss_bc03_color,100)
       oplot, bc03_supersolar.zobj, bc03_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(bc03_psym,thick=6), color=fsc_color(bc03_color,100)
    endelse

; basti-ss
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ss_supersolar.zobj, sdss_basti_ss_supersolar.mass_50, sdss_basti_ss_supersolar.mass_err, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100), $
         errcolor=fsc_color(sdss_basti_ss_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ss_supersolar.zobj, basti_ss_supersolar.mass_50, basti_ss_supersolar.mass_err, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100), $
         errcolor=fsc_color(basti_ss_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ss_supersolar.zobj, sdss_basti_ss_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ss_psym,thick=6), color=fsc_color(sdss_basti_ss_color,100)
       oplot, basti_ss_supersolar.zobj, basti_ss_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(basti_ss_psym,thick=6), color=fsc_color(basti_ss_color,100)
    endelse

; basti-ae
    if keyword_set(showerr) then begin
       oploterror, sdss_basti_ae_supersolar.zobj, sdss_basti_ae_supersolar.mass_50, sdss_basti_ae_supersolar.mass_err, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100), $
         errcolor=fsc_color(sdss_basti_ae_color,100), symsize=symsize1, errthick=4
       oploterror, basti_ae_supersolar.zobj, basti_ae_supersolar.mass_50, basti_ae_supersolar.mass_err, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100), $
         errcolor=fsc_color(basti_ae_color,100), symsize=symsize1, errthick=4
    endif else begin
       oplot, sdss_basti_ae_supersolar.zobj, sdss_basti_ae_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(sdss_basti_ae_psym,thick=6), color=fsc_color(sdss_basti_ae_color,100)
       oplot, basti_ae_supersolar.zobj, basti_ae_supersolar.mass_50, symsize=symsize1, $
         psym=symcat(basti_ae_psym,thick=6), color=fsc_color(basti_ae_color,100)
    endelse

    im_legend, ['BC03','BaSTI-ss','BaSTI-ae'], /right, /top, box=0, $
      color=[bc03_color,basti_ss_color,basti_ae_color], $
      psym=[bc03_psym,basti_ss_psym,basti_ae_psym], $
      symthick=4, charsize=1.5
    xyouts, 0.8, 12.2, 'Z = 0.05', align=0.5, /data, charsize=1.5

    im_plotconfig, psfile=psfile, /psclose, /gzip

stop    
    
return
end
    

; read the 37 kpc K-band photometry and compute the aperture
; correction
    sample = rsex('bcgs_sample_v3.sex')
    phot = mrdfits('bcgs_photometry_v3.fits.gz',1)
    all = djs_readlines('all.Ks.out')
    ngal = n_elements(all)
    ktot = dblarr(ngal)
;   for ii = 0, ngal-1 do ktot[ii] = (strsplit(all[ii],' ',/extract))[3]
    for ii = 0, ngal-1 do ktot[ii] = (strsplit(all[ii],' ',/extract))[38]
    niceprint, lindgen(n_elements(sample)), sample.z, sample.dec, ktot

    im_plothist, 10^(-0.4*(phot.ks_mag_aper_08-ktot))
    
