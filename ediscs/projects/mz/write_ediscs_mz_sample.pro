pro write_ediscs_mz_sample, ancillary, ispec, write=write
; jm08apr07nyu - written, based on WRITE_AGES_MZ_SAMPLE

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()

    lsun = 3.826D33
    hahb = 2.86 ; case B
    haoii = 1.0 ; assumed intrinsic ratio
    haconst = 7.9D-42 ; K98 conversion L(Ha) --> SFR(Ha) (Salpeter 0.1-100)
    absmsun_b = k_solar_magnitudes(filterlist='bessell_B.par')

    mzpath = ediscs_path(/projects)+'mz/'
    ages_mzpath = ages_path(/projects)+'mz/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = mzpath+'write_ediscs_mz_sample.log'
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

    plotsym, 0, 0.1, /fill
    charsize_6 = 1.6

    ewhbcut = 1.0
    ewoiicut = 6.0 ; 3.0
    snrcut = 5.0
    fwhm2sig = 2.35
    broad_sigmacut = 700.0
    
    snrrange = [0.01,200.0]
    ewrange = [0.001,500.0]
    ewrange2 = [0.5,500.0]
    kmsrange = [0.1,2000.0]

; read the SDSS coefficients for computing "alpha"    

    if (file_test(ages_mzpath+'sdss_alphafit.fits') eq 0L) then begin
       splog, '"Alpha" file '+ages_mzpath+'sdss_alphafit.fits not found.'
       return
    endif

    sdss_alphafit = mrdfits(ages_mzpath+'sdss_alphafit.fits',1,/silent)
    alpha_d4000_coeff = [sdss_alphafit[0].alpha_coeff]
    alpha_gr_coeff = [sdss_alphafit[1].alpha_coeff]
    
; ---------------------------------------------------------------------------    
; read the data
; ---------------------------------------------------------------------------    

    splog, 'Reading the data...'
    if (n_elements(ancillary) eq 0L) then ancillary = read_ediscs(/ancillary)
    if (n_elements(ispec) eq 0L) then begin
       ispec1 = read_ediscs(/specfit)
       ispecmore = {$
         final_class:         '?', $
         r23branch:           '?', $
         ewalpha_gr:       -999.0, $
         ewalpha_d4000:    -999.0, $
         oii_lum:          -999.0, $
;        oii_lum_err:      -999.0, $
         hb_lum:           -999.0, $
;        hb_lum_err:       -999.0, $
         sfr_oii:          -999.0, $
;        sfr_oii_err:      -999.0, $
         sfr_hb:           -999.0, $
;        sfr_hb_err:       -999.0, $
         sfr_oii_mk06:     -999.0, $
;        sfr_oii_mk06_err: -999.0, $
         sfr_hb_mk06:      -999.0}
;        sfr_hb_mk06_err:  -999.0}
       ispecmore = replicate(ispecmore,n_elements(ispec1))
       ispec = struct_addtags(ispec1,ispecmore)
    endif
    ngalaxy = n_elements(ancillary)

; ---------------------------------------------------------------------------    
; parent sample cuts; 0.033<z<0.25; 14.5<r<17.77; well measured g,r,i
; magnitudes (for k-corrections and masses); the cuts on the
; photometry, INFIBER, and MASS are *very* soft, removing <0.05% of
; the sample; note that we cut on apparent magnitudes to ensure that
; the k-corrections are reliable
; ---------------------------------------------------------------------------    

    survey_area = 1.0*!dtor^2.0 ; 3.05D-4 sr/deg^2 [sr]
    vname = 'default.nolines'
    survey_filter = 'FORS_I_ccd_atm.par'
    sample_zmin = 0.01
    sample_zmax = 1.0
    ediscs_mlimit = [15.0,24.0]
    
    parent = where((ancillary.z gt sample_zmin) and (ancillary.z lt sample_zmax),nparent)
    splog, 'Parent sample: '+string(nparent,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      strtrim(string(100.0*nparent/float(ngalaxy),format='(F12.1)'),2)+'%).'

    parent_ancillary = ancillary[parent]
    parent_ispec = ispec[parent]

; statistics    
    
    zstats = im_stats(parent_ancillary.z)
    splog, '   Redshift      : ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; ---------------------------------------------------------------------------
; define the mass-metallicity sample; the cuts on positive EW are
; *very* soft: the more generous cut is that we include objects that
; are within 1-sigma of the EWCUT; the upper cut on EW([OII]) removes
; one object with a crazy measurement; the cuts on the error being
; positive removes just a handful of low-redshift objects where [OII]
; is very close to the blue wavelength cut-off; in EDISCS_MZ_LOG12OH we
; compute the metallicities using the equivalent widths
; ---------------------------------------------------------------------------

    mz1 = where($
      (parent_ispec.h_beta_ew[0] gt 0.0) and (parent_ispec.oii_3727_ew[0] gt 0.0) and $
      (parent_ispec.h_beta_ew[0]+parent_ispec.h_beta_ew[1] gt ewhbcut) and $
      (parent_ispec.oii_3727_ew[0]+parent_ispec.oii_3727_ew[1] gt ewoiicut) and $
      (parent_ispec.h_beta_ew[1] gt 0.0) and (parent_ispec.oii_3727_ew[1] gt 0.0) and $
      (parent_ispec.oiii_5007_ew[1] gt 0.0),nmz1)

    splog, 'MZ sample: '+string(nmz1,format='(I0)')+'/'+string(nparent,format='(I0)')+' ('+$
      strtrim(string(100.0*nmz1/nparent,format='(F12.1)'),2)+'%).'

    mz_ancillary = parent_ancillary[mz1]
    mz_ispec = parent_ispec[mz1]
    mz_ancillary = parent_ancillary[mz1]
    mz_ispec = parent_ispec[mz1]

; statistics    
    
    zstats = im_stats(mz_ancillary.z)
    splog, '   Redshift      : ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; compare the two MZ samples

;   im_plothist, mz_ancillary.ugriz_absmag[2], bin=0.1, charsize=1.5
;   im_plothist, ancillary[parent[mz2]].ugriz_absmag[2], bin=0.1, /overplot, line=2, thick=2
;
;   im_plothist, mz_ancillary.ugriz_absmag[1]-mz_ancillary.ugriz_absmag[2], $
;     bin=0.05, charsize=1.5, xr=[-0.2,1.3]
;   im_plothist, ancillary[parent[mz2]].ugriz_absmag[1]-ancillary[parent[mz2]].ugriz_absmag[2], $
;     bin=0.05, /overplot, line=2, thick=2

; ---------------------------------------------------------------------------
; compute emission-line star-formation rates
; ---------------------------------------------------------------------------

    splog, 'Computing emission-line star-formation rates.'

    loglb = alog10(10.0^(-0.4*(mz_ancillary.ubvrijhk_absmag[1]-absmsun_b)))
    
    dlum = dluminosity(mz_ancillary.z,/cm)
    hb_lum_factor = mz_ancillary.cflux_4861*(4.0*!dpi*dlum^2.0)/lsun
    oii_lum_factor = mz_ancillary.cflux_3727*(4.0*!dpi*dlum^2.0)/lsun
    
    mz_ispec.hb_lum      = alog10(mz_ispec.h_beta_ew[0]*hb_lum_factor)
    mz_ispec.oii_lum     = alog10(mz_ispec.oii_3727_ew[0]*oii_lum_factor)
;   mz_ispec.hb_lum_err  = mz_ispec.h_beta_ew[1]/mz_ispec.h_beta_ew[0]/alog(10.0)
;   mz_ispec.oii_lum_err = mz_ispec.oii_3727_ew[1]/mz_ispec.oii_3727_ew[0]/alog(10.0)

    mz_ispec.sfr_hb = mz_ispec.hb_lum+alog10(lsun)+alog10(hahb)+alog10(haconst)
    mz_ispec.sfr_oii = mz_ispec.oii_lum+alog10(lsun)+alog10(haoii)+alog10(haconst)

    mz_ispec.sfr_hb_mk06 = hb_sfr(loglb,mz_ispec.hb_lum+alog10(lsun),sfr_err=sfr_hb_err,/log)
    mz_ispec.sfr_oii_mk06 = oii_sfr(loglb,mz_ispec.oii_lum+alog10(lsun),sfr_err=sfr_oii_err,/log)

; ---------------------------------------------------------------------------
; assign EWALPHA values
; ---------------------------------------------------------------------------

    splog, 'Assigning EW alpha values.'

; note that the alpha-Dn(4000) relation is only strictly calibrated
; for 1.1<Dn(4000)<1.8, and the alpha-^{0.1}(g-r) relation is only
; calibrated for 0.15<^{0.1}(g-r)<1.0
    
    d4000 = mz_ispec.d4000_narrow_model[0]
    mz_ispec.ewalpha_d4000 = poly(d4000,alpha_d4000_coeff)

    gr = mz_ancillary.ugriz_absmag[1]-mz_ancillary.ugriz_absmag[2]
    mz_ispec.ewalpha_gr = poly(gr,alpha_gr_coeff)

;   plot, mz_ispec.ewalpha_d4000, mz_ispec.ewalpha_gr, ps=3, $
;     xr=[0.7,1.5], yr=[0.7,1.5] , xsty=3, ysty=3, charsize=2
;   oplot, !x.crange, !y.crange

; ---------------------------------------------------------------------------
; assign R23 branches
; ---------------------------------------------------------------------------

    splog, 'Assigning all objects to the upper R23 branch.'
    mz_ispec.r23branch = 'U'

;;    splog, 'Assigning R23 branches based on [N II]/Ha.'
;;    r23branch_niiha = where((mz_ispec.h_alpha[0] gt 0.0) and (mz_ispec.h_alpha[1] gt 0.0) and $
;;      (mz_ispec.nii_6584[0] gt 0.0) and (mz_ispec.nii_6584[1] gt 0.0) and $
;;      (mz_ispec.h_alpha[0]/mz_ispec.h_alpha[1] gt snrcut),nr23branch_niiha,$
;;      comp=r23branch_default,ncomp=nr23branch_default)
;;    if (nr23branch_default ne 0L) then mz_ispec[r23branch_default].r23branch = 'U'
;;
;;;   niioii = mz_ispec[r23branch_niiha].nii_6584[0]/mz_ispec[r23branch_niiha].oii_3727[0]
;;    niiha = mz_ispec[r23branch_niiha].nii_6584[0]/mz_ispec[r23branch_niiha].h_alpha[0]
;;    niiha_lower = where(niiha lt 10^(-1.0),nniiha_lower,comp=niiha_upper,ncomp=nniiha_upper)
;;    if (nniiha_lower ne 0L) then mz_ispec[r23branch_niiha[niiha_lower]].r23branch = 'L'
;;    if (nniiha_upper ne 0L) then mz_ispec[r23branch_niiha[niiha_upper]].r23branch = 'U'
;;
;;    splog, '   No assignment possible (default "U" assumed): '+string(nmz1-nr23branch_niiha,format='(I0)')+'/'+$
;;      string(nmz1,format='(I0)')+' ('+strtrim(string(100.0*(nmz1-nr23branch_niiha)/nmz1,format='(F12.1)'),2)+'%).'
;;    splog, '   Lower branch: '+string(nniiha_lower,format='(I0)')+'/'+string(nr23branch_niiha,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nniiha_lower/nr23branch_niiha,format='(F12.1)'),2)+'%).'
;;    splog, '   Upper branch: '+string(nniiha_upper,format='(I0)')+'/'+string(nr23branch_niiha,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nniiha_upper/nr23branch_niiha,format='(F12.1)'),2)+'%).'
;;    print

; ---------------------------------------------------------------------------
; identify AGN
; ---------------------------------------------------------------------------
;;
;;; BPT diagram    
;;
;;    splog, 'Classifying...' & print
;;    iclass = iclassification(mz_ispec,ratios=iratios,/kauffman,$
;;      snrcut_class=snrcut,/silent,/doplot)
;;
;;    canclass = where((strtrim(iclass.bpt_pure_nii_class,2) ne 'Unknown'),$
;;      ncanclass,comp=noclass,ncomp=nnoclass)
;;    if (nnoclass ne 0L) then mz_ispec[noclass].final_class = 'BPT_NOCLASS'
;;    agn = where(strtrim(iclass.bpt_pure_nii_class,2) eq 'AGN',nagn)
;;    if (nagn ne 0L) then mz_ispec[agn].final_class = 'BPT_AGN'
;;    hii = where(strtrim(iclass.bpt_pure_nii_class,2) eq 'HII',nhii)
;;    if (nhii ne 0L) then mz_ispec[hii].final_class = 'BPT_HII'    
;;    hiiplus = where((strtrim(iclass.bpt_pure_nii_class,2) eq 'HII') or $
;;      (strtrim(iclass.bpt_pure_nii_class,2) eq 'Unknown'),nhiiplus)
;;    
;;    splog, 'BPT/Unclassified: '+string(nnoclass,format='(I0)')+'/'+string(nmz1,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nnoclass/nmz1,format='(F12.1)'),2)+'%).'
;;    splog, 'BPT/Classified  : '+string(ncanclass,format='(I0)')+'/'+string(nmz1,format='(I0)')+' ('+$
;;      strtrim(string(100.0*ncanclass/nmz1,format='(F12.1)'),2)+'%).'
;;    splog, '   BPT/HII: '+string(nhii,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nhii/ncanclass,format='(F12.1)'),2)+'%).'
;;    splog, '   BPT/AGN: '+string(nagn,format='(I0)')+'/'+string(ncanclass,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nagn/ncanclass,format='(F12.1)'),2)+'%).'
;;    splog, 'BPT/HII+Unclassified: '+string(nhiiplus,format='(I0)')+'/'+string(nmz1,format='(I0)')+' ('+$
;;      strtrim(string(100.0*nhiiplus/nmz1,format='(F12.1)'),2)+'%).'
;;    print
;;
;;    mzhii_ancillary = mz_ancillary[hii]
;;    mzhiiplus_ancillary = mz_ancillary[hiiplus]
;;    mzagn_ancillary = mz_ancillary[agn]
;;    mzhii_ispec = mz_ispec[hii]
;;    mzhiiplus_ispec = mz_ispec[hiiplus]
;;    mzagn_ispec = mz_ispec[agn]

; ---------------------------------------------------------------------------    
; compute the dust attenuation
; ---------------------------------------------------------------------------    

;   splog, 'Computing the dust attenuation...'
;   mzhii_ispec_nodust = iunred_linedust(mzhii_ispec,snrcut=snrcut,/silent)
;   mzhiiplus_ispec_nodust = iunred_linedust(mzhiiplus_ispec,snrcut=snrcut,/silent)

; ###########################################################################
; interlude: sample selection plots; put them in a single PS file
; ###########################################################################

    ewhb = parent_ispec.h_beta_ew[0]
    ewhb_err = parent_ispec.h_beta_ew[1]
    snrhb = parent_ispec.h_beta_ew[0]/parent_ispec.h_beta_ew[1]

    ewoii = parent_ispec.oii_3727_ew[0]
    ewoii_err = parent_ispec.oii_3727_ew[1]
    snroii = parent_ispec.oii_3727_ew[0]/parent_ispec.oii_3727_ew[1]

    ewoiii = parent_ispec.oiii_5007_ew[0]
    ewoiii_err = parent_ispec.oiii_5007_ew[1]
    snroiii = parent_ispec.oiii_5007_ew[0]/parent_ispec.oiii_5007_ew[1]

    ewaxis = findgen((alog10(1E3)-alog10(1E-3))/0.05)*0.05+alog10(1E-3)

; Hb    
    
    hb_good = where((ewhb gt 0.0) and (snrhb gt 0.0))
    hb_coeff = robust_linefit(alog10(ewhb[hb_good]),alog10(snrhb[hb_good]),/bisect,hb_yfit,hb_sig)
    hbaxis = poly(ewaxis,hb_coeff)

    ewhbcut_1 = interpol(10.0^ewaxis,10.0^hbaxis,1.0)
    ewhbcut_3 = interpol(10.0^ewaxis,10.0^hbaxis,3.0)
    ewhbcut_5 = interpol(10.0^ewaxis,10.0^hbaxis,5.0)
    splog, 'EW(Hb)      at S/N EW(Hb)      = 1, 3, 5: ', ewhbcut_1, ewhbcut_3, ewhbcut_5

; [O II]
    
    oii_good = where((ewoii gt 0.0) and (snroii gt 0.0))
    oii_coeff = robust_linefit(alog10(ewoii[oii_good]),alog10(snroii[oii_good]),/bisect,oii_yfit,oii_sig)
    oiiaxis = poly(ewaxis,oii_coeff)

    ewoiicut_1 = interpol(10.0^ewaxis,10.0^oiiaxis,1.0)
    ewoiicut_3 = interpol(10.0^ewaxis,10.0^oiiaxis,3.0)
    ewoiicut_5 = interpol(10.0^ewaxis,10.0^oiiaxis,5.0)
    splog, 'EW([O II])  at S/N EW([O II])  = 1, 3, 5: ', ewoiicut_1, ewoiicut_3, ewoiicut_5
    
; [O III]

    oiii_good = where((ewoiii gt 0.0) and (snroiii gt 0.0))
    oiii_coeff = robust_linefit(alog10(ewoiii[oiii_good]),alog10(snroiii[oiii_good]),/bisect,oiii_yfit,oiii_sig)
    oiiiaxis = poly(ewaxis,oiii_coeff)

    ewoiiicut_1 = interpol(10.0^ewaxis,10.0^oiiiaxis,1.0)
    ewoiiicut_3 = interpol(10.0^ewaxis,10.0^oiiiaxis,3.0)
    ewoiiicut_5 = interpol(10.0^ewaxis,10.0^oiiiaxis,5.0)
    splog, 'EW([O III]) at S/N EW([O III]) = 1, 3, 5: ', ewoiiicut_1, ewoiiicut_3, ewoiiicut_5

; index for testing    

    hb_mz_indx = where((ewhb+ewhb_err gt ewhbcut) and (ewhb/ewhb_err gt 1.0))
    mz_indx = where((ewhb+ewhb_err gt ewhbcut) and (ewhb/ewhb_err gt 1.0) and $
      (ewoii+ewoii_err gt ewoiicut) and (ewoii/ewoii_err gt 1.0))

; now make the plot    
    
    if keyword_set(write) then begin
       psname = mzpath+'ediscs_mz_sample_selection.ps'
       dfpsplot, psname, /color, /square
    endif

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.6,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.1, $
      position=pos, /normal

; ------------------------------------------------------------
; EW(Hb) vs S/N EW(Hb)
 
    xtitle = 'EW(H\beta) [\AA]'
    ytitle = 'S/N EW(H\beta)'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
    djs_oplot, ewhb, snrhb, psym=3, color='grey'
    djs_oplot, 10.0^ewaxis, 10.0^hbaxis, line=0, thick=postthick2

    djs_oplot, ewhbcut_1*[1,1], [10^!y.crange[0],1.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_3*[1,1], [10^!y.crange[0],3.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_5*[1,1], [10^!y.crange[0],5.0], line=2, thick=postthick2

    djs_oplot, [10^!x.crange[0],ewhbcut_1], 1.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_3], 3.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_5], 5.0*[1,1], line=2, thick=postthick2

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; EW([O II]) vs S/N EW([OII])

    xtitle = 'EW([O II]) [\AA]'
    ytitle = 'S/N EW([O II])'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
    djs_oplot, ewoii, snroii, psym=3, color='grey'
    djs_oplot, ewoii[hb_mz_indx], snroii[hb_mz_indx], psym=3, color='dark green'
    djs_oplot, 10.0^ewaxis, 10.0^oiiaxis, line=0, thick=postthick2
    legend, textoidl('EW(H\beta) > '+string(ewhbcut,format='(F3.1)')), /left, /top, symsize=3.0, $
      box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

    djs_oplot, ewoiicut_1*[1,1], [10^!y.crange[0],1.0], line=2, thick=postthick2
    djs_oplot, ewoiicut_3*[1,1], [10^!y.crange[0],3.0], line=2, thick=postthick2
    djs_oplot, ewoiicut_5*[1,1], [10^!y.crange[0],5.0], line=2, thick=postthick2

    djs_oplot, [10^!x.crange[0],ewoiicut_1], 1.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewoiicut_3], 3.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewoiicut_5], 5.0*[1,1], line=2, thick=postthick2

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; EW(Hb) vs S/N EW([OII])

    xtitle = 'EW(H\beta) [\AA]'
    ytitle = 'S/N EW([O II])'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
    djs_oplot, ewhb, snroii, psym=3, color='grey'
    djs_oplot, ewhb[hb_mz_indx], snroii[hb_mz_indx], psym=3, color='dark green'
    legend, textoidl('EW(H\beta) > '+string(ewhbcut,format='(F3.1)')), /left, /top, symsize=3.0, $
      box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

    djs_oplot, ewhbcut_1*[1,1], [10^!y.crange[0],1.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_3*[1,1], [10^!y.crange[0],3.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_5*[1,1], [10^!y.crange[0],5.0], line=2, thick=postthick2

    djs_oplot, [10^!x.crange[0],ewhbcut_1], 1.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_3], 3.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_5], 5.0*[1,1], line=2, thick=postthick2

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; EW([O III]) vs S/N EW([OIII])

    xtitle = 'EW([O III]) [\AA]'
    ytitle = 'S/N EW([O III])'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
    djs_oplot, ewoiii, snroiii, psym=3, color='grey'
    djs_oplot, ewoiii[mz_indx], snroiii[mz_indx], psym=3, color='dark green'
    djs_oplot, 10.0^ewaxis, 10.0^oiiiaxis, line=0, thick=postthick2
    legend, textoidl('EW(H\beta) > '+string(ewhbcut,format='(F3.1)')+'; '+$
      'EW([O II]) > '+string(ewoiicut,format='(F3.1)')), /left, /top, symsize=3.0, $
      box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

    djs_oplot, ewoiiicut_1*[1,1], [10^!y.crange[0],1.0], line=2, thick=postthick2
    djs_oplot, ewoiiicut_3*[1,1], [10^!y.crange[0],3.0], line=2, thick=postthick2
    djs_oplot, ewoiiicut_5*[1,1], [10^!y.crange[0],5.0], line=2, thick=postthick2

    djs_oplot, [10^!x.crange[0],ewoiiicut_1], 1.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewoiiicut_3], 3.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewoiiicut_5], 5.0*[1,1], line=2, thick=postthick2

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

; ------------------------------------------------------------
; EW(Hb) vs S/N EW([OIII])

    xtitle = 'EW(H\beta) [\AA]'
    ytitle = 'S/N EW([O III])'
    xrange = ewrange
    yrange = snrrange

    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
    djs_oplot, ewhb, snroiii, psym=3, color='grey'
    djs_oplot, ewhb[mz_indx], snroiii[mz_indx], psym=3, color='dark green'
    legend, textoidl('EW(H\beta) > '+string(ewhbcut,format='(F3.1)')+'; '+$
      'EW([O II]) > '+string(ewoiicut,format='(F3.1)')), /left, /top, symsize=3.0, $
      box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')

    djs_oplot, ewhbcut_1*[1,1], [10^!y.crange[0],1.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_3*[1,1], [10^!y.crange[0],3.0], line=2, thick=postthick2
    djs_oplot, ewhbcut_5*[1,1], [10^!y.crange[0],5.0], line=2, thick=postthick2

    djs_oplot, [10^!x.crange[0],ewhbcut_1], 1.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_3], 3.0*[1,1], line=2, thick=postthick2
    djs_oplot, [10^!x.crange[0],ewhbcut_5], 5.0*[1,1], line=2, thick=postthick2

    if (not keyword_set(write)) then begin
       splog, 'Press any key to continue.'
       cc = get_kbrd(1)
    endif

;; ------------------------------------------------------------
;; EW(Hb) vs error in 12+log(O/H)
;    
;    xtitle = 'EW(H\beta)'
;    ytitle = '\sigma[12+log(O/H)]'
;    xrange = [0.01,500]
;    yrange = [0.001,10.0]
; 
;    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, xrange=xrange, yrange=yrange, $
;      xtitle=xtitle, ytitle=ytitle, charsize=charsize_6, charthick=postthick1, $
;      xthick=postthick1, ythick=postthick1, position=pos[*,0], /xlog, /ylog
;    djs_oplot, ewhb[ohgood], oherr[ohgood], psym=3, color='grey'
;    djs_oplot, ewhb[mz_indx_ohgood], oherr[mz_indx_ohgood], psym=3, color='dark green'
;     legend, textoidl('EW(H\beta) > '+string(ewhbcut,format='(F3.1)')+'; '+$
;       'EW([O II]) > '+string(ewoiicut,format='(F3.1)')), /left, /top, symsize=3.0, $
;       box=0, charsize=1.5, charthick=postthick1, psym=3, color=djs_icolor('dark green')
;    
;    if (not keyword_set(write)) then begin
;       splog, 'Press any key to continue.'
;       cc = get_kbrd(1)
;    endif
 
; ------------------------------------------------------------

    if keyword_set(write) then begin
       dfpsclose
       spawn, 'gzip -f '+psname
    endif

; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin

; ancillary data
       
       parentfile_ancillary = mzpath+'ediscs_parent_ancillary.fits'
       splog, 'Writing '+parentfile_ancillary
       mwrfits, parent_ancillary, parentfile_ancillary, /create

       mzfile_ancillary = mzpath+'ediscs_mz_ancillary.fits'
       splog, 'Writing '+mzfile_ancillary
       mwrfits, mz_ancillary, mzfile_ancillary, /create

; emission-line data       

       mzfile_ispec = mzpath+'ediscs_mz_ispec.fits'
       splog, 'Writing '+mzfile_ispec
       mwrfits, mz_ispec, mzfile_ispec, /create

       splog, 'Total time = '+strtrim(string((systime(1)-stime0)/60.0,format='(F12.1)'),2)+' minutes.'
       splog, /close

; finally compute and write out the oxygen abundances

       splog, 'Computing oxygen abundances.'
       mz_log12oh, mz_ispec, nmonte=500L, snrcut=0.0, $
         logfile=mzpath+'mz_log12oh.ediscs.log', $
         outfile_root=mzpath+'ediscs_mz_log12oh', /write
       
    endif

stop

return
end
