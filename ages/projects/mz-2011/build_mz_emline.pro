function doclass, mz_ispec, iratios=iratios, snrcut1=snrcut1, doplot=doplot
; classify using various techniques

    splog, 'Classifying...'
    ngal = n_elements(mz_ispec)

; BPT diagram    
    class = replicate({agn: 0, bpt_d: 0.0, bpt_phi: 0.0, class: ''},ngal)
    iclass = iclassification(mz_ispec,ratios=iratios,$
      snrcut_class=snrcut1,silent=0,doplot=doplot)
    class.bpt_d = iratios.bpt_d
    class.bpt_phi = iratios.bpt_phi
    class.class = strtrim(iratios.final_class,2)

; rename the classes to make the matching easier    
    agn = where(strtrim(class.class,2) eq 'AGN')
    if (agn[0] ne -1) then class[agn].class = 'BPT_AGN'
    sfagn = where(strtrim(class.class,2) eq 'SF/AGN')
    if (sfagn[0] ne -1) then class[sfagn].class = 'BPT_AGN'
    sf = where(strtrim(class.class,2) eq 'SF')
    if (sf[0] ne -1) then class[sf].class = 'BPT_SF'
;   unknown = where(strtrim(class.class,2) eq 'UNKNOWN')
;   if (unknown[0] ne -1) then class[unknown].class = 'UNKNOWN'
    
return, class
end

pro build_mz_emline, sdss=sdss, clobber=clobber, doplot=doplot
; jm09mar18nyu - build the AGES and SDSS emission-line samples; also
; classify into SF, AGN, SF/AGN, and UNKNOWN

    mzpath = mz_path()
    catpath = ages_path(/catalogs)

    parent = read_mz_sample(sdss=sdss,/parent)
    mass = read_mz_sample(sdss=sdss,/mass)

; see MZPLOT_SAMPLE for details    
    ages_hbcut = mz_hbcut()
    ages_ewcut = 0.0
    ages_snrcut = 1.0
    ages_oiihbcut = -0.3

    sdss_hbcut = mz_hbcut(/sdss)
    sdss_ewcut = 0.0
    sdss_snrcut = 1.0
    sdss_oiihbcut = -0.3

    if keyword_set(sdss) then begin ; SDSS
       splog, file=mzpath+'build_mz_emline.sdss.log'
       splog, '#########################'
       splog, 'Building SDSS emission-line sample'

; POSTSTR/35: -17>Mr>-24; 0.033<z<0.25; there are a handful (9) of
; objects with *ridiculous [OII] EWs (>10,000 A) - toss them here!
       vagc = mz_get_vagc(sample=sample,letter=letter,poststr=poststr)
       ispec = read_vagc_garching(sample=sample,$
         letter=letter,poststr=poststr,/ispec)
       ispec = ispec[parent.vagc_object_position]
       
       mz = where($
         (ispec.h_beta[0] gt sdss_hbcut) and $
         (ispec.h_beta[0]/ispec.h_beta[1] gt sdss_snrcut) and $ ; nominal cut
;        (ispec.h_beta_ew[0] gt sdss_ewcut) and $
;        (ispec.h_beta_sigma[0] lt sdss_sigmacut) and $
         (ispec.oii_3727_ew[0] lt 1000.0) and $
         (ispec.oii_3727[1] ne -2.0) and $
         (ispec.oiii_5007[1] ne -2.0) and $ ; nominal
         (ispec.nii_6584[1] ne -2.0) and $  ; nominal
         (ispec.h_alpha[1] ne -2.0) and $   ; nominal
         (ispec.oii_3727[0] gt 10^sdss_oiihbcut*ispec.h_beta[0]),ngal)
       splog, 'Emission-line sample = ', ngal

       mz_ispec = ispec[mz]
       mz_parent = parent[mz]
       mz_mass = mass[mz]

       splog, '[OIII] upper limits = ', total(mz_ispec.oiii_5007[1] eq -1), $
         total(mz_ispec.oiii_5007[1] eq -1)/float(ngal)

; classify into SF and AGN using various techniques

; 1) BPT diagram, including upper limits; do not use the pure [NII]/Ha
; classification since we deal with upper limits below
       iclass = iclassification(mz_ispec,ratios=iratios1,$
         snrcut_class=0.0,silent=0,doplot=doplot,/noniihaclass)
;      class = doclass(mz_ispec,iratios=iratios,$
;        snrcut1=sdss_snrcut,doplot=doplot)
       bpt_agn = (strtrim(iratios1.final_class,2) eq 'AGN') or $
         (strtrim(iratios1.final_class,2) eq 'SF/AGN')
       bpt_sf = strtrim(iratios1.final_class,2) eq 'SF'
       
       iratios = mzbpt_upper_limits(iratios1,mz_ispec,doplot=doplot)
       bpt_sf_limit = strmatch(iratios.final_class,'BPT_SF_*LIMIT*')
       bpt_agn_limit = strmatch(iratios.final_class,'BPT_AGN_*LIMIT*')

       splog, 'BPT-SF ', total(bpt_sf), total(bpt_sf)/float(ngal)
       splog, 'BPT-AGN ', total(bpt_agn), total(bpt_agn)/float(ngal)
       splog, 'BPT-SF-LIMITS ', total(bpt_sf_limit), total(bpt_sf_limit)/float(ngal)
       splog, 'BPT-AGN-LIMITS ', total(bpt_agn_limit), total(bpt_agn_limit)/float(ngal)

; note: a small fraction of objects (2080/198024=1.05%) cannot be
; classified using the above techniques either because their position
; in the BPT diagram was ambiguous or because H-alpha was not detected
; (1289/2080=62%)
       nobptclass = (strtrim(iratios.final_class,2) eq 'UNKNOWN')
       noha = (mz_ispec[where(nobptclass)].h_alpha[1] ne -2.0) and $
         (iratios[where(nobptclass)].nii_ha lt -900) and $
         (iratios[where(nobptclass)].nii_ha_limit lt -900)

       splog, 'No BPT classification = ', total(nobptclass), total(nobptclass)/float(ngal)
       splog, '  Of these, no H-alpha = ', total(noha), total(noha)/total(nobptclass)
       
; Yan diagram; all the [OIII]/Hb upper limits are consistent with
; being star-forming, while the BPT-unclassified objects are
; essentially all Yan-AGN
       ub = mz_parent.k_ubvrijhk_absmag_00[0]-mz_parent.k_ubvrijhk_absmag_00[1]
       yanagn = (iratios.oiii_hb ge ((1.4-1.2*ub)>(-0.1))) or $
         (iratios.oiii_hb_limit ge ((1.4-1.2*ub)>(-0.1)))
       yanandnobpt = (yanagn eq 1) and (nobptclass)
       yanandoiiilimit = (yanagn eq 1) and (iratios.oiii_hb_limit gt -900)

       splog, 'Yan AGN', total(yanagn), total(yanagn)/float(ngal)
       splog, '  Of these, no BPT class ', total(yanandnobpt), total(yanandnobpt)/float(ngal)
       splog, '  Of these, [OIII]/Hb limit ', total(yanandoiiilimit), total(yanandoiiilimit)/float(ngal)
       
; debugging plot
       if keyword_set(doplot) then begin
          det = where(iratios.oiii_hb gt -900.0,comp=lim)
          djs_plot, ub[det], iratios[det].oiii_hb, psym=6, xr=[0,2], sym=0.1, yr=[-2,2]
          djs_oplot, ub[lim], iratios[lim].oiii_hb_limit, psym=6, sym=0.1, color='yellow'
          ubaxis = range(0.0,2.0,500)
          djs_oplot, ubaxis, (1.4-1.2*ubaxis)>(-0.1), line=0, thick=7
          wdet = where(yanagn[det] eq 1)
          wlim = where(yanagn[lim] eq 1)
          djs_oplot, ub[det[wdet]], iratios[det[wdet]].oiii_hb, psym=6, sym=0.1, color='cyan'
;         djs_oplot, ub[wlim], iratios[wlim].oiii_hb_limit, psym=6, sym=0.2, color='yellow'
       endif
       
; pack it all into a structure and then write out
       class = replicate({agn: -1, bpt_agn: -1, yan_agn: 0},ngal)
       class = struct_addtags(class,im_struct_trimtags(iratios,$
         select=['bpt_d','bpt_phi','nii_ha','nii_ha_err','oiii_hb',$
         'oiii_hb_err','nii_ha_limit','oiii_hb_limit'],$
         newtags=['bpt_d','bpt_phi','bpt_nii_ha','bpt_nii_ha_err','bpt_oiii_hb',$
         'bpt_oiii_hb_err','bpt_nii_ha_limit','bpt_oiii_hb_limit']))

       class[where(bpt_agn or bpt_agn_limit)].bpt_agn = 1 ; 1=agn, 0=sf, -1=unclassified
       class[where(bpt_sf or bpt_sf_limit)].bpt_agn = 0
       class.yan_agn = yanagn

; final classifications: final class is AGN if it's been classified as
; an AGN using any of the other criteria above, unless it's been
; classified as star-forming using the BPT diagram
       class.agn = (class.bpt_agn eq 1) or (class.yan_agn eq 1)
       reset = where(class.agn eq 1 and class.bpt_agn eq 0,nreset)
       if (nreset ne 0L) then class[reset].agn = 0

; some statistics
;      nclass = total(class.bpt_agn eq -1)
;      yclass = total(class.bpt_agn ne -1)
;      yagn = total((class.bpt_agn ne -1) and (class.bpt_agn eq 1))
;      ysf = total((class.bpt_agn ne -1) and (class.bpt_agn eq 0))
       splog, 'Final statistics:'
       splog, 'AGN: ', total(class.agn eq 1), ngal, total(class.agn eq 1)/float(ngal)
       splog, 'SF : ', total(class.agn eq 0), ngal, total(class.agn eq 0)/float(ngal)
       
; store the final classifications and then write out
       mz_ispec = struct_addtags(temporary(mz_ispec),class)
       hii = where(mz_ispec.agn eq 0,nhii,comp=agn,ncomp=nagn)
       
;      djs_plot, class[agn].bpt_d, class[agn].bpt_phi, psym=3, xr=[-0.1,2.5], $
;        yr=[-90,90], xsty=1, ysty=1, color='orange'
;      djs_oplot, class[hii].bpt_d, class[hii].bpt_phi, psym=3, color='cyan'
;      djs_oplot, !x.crange, 25*[1,1]
;
;      ww = where(iratios[agn].bpt_phi lt -30)
;      djs_plot, iratios.nii_ha, iratios.oiii_hb, psym=3, xr=[-1.5,1]                 
;      djs_oplot, iratios[agn[ww]].nii_ha, iratios[agn[ww]].oiii_hb, psym=3, color='red'
       
       im_mwrfits, mz_parent, mzpath+'sdss_mz_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec, mzpath+'sdss_mz_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass, mzpath+'sdss_mz_mass.fits', clobber=clobber

       im_mwrfits, mz_parent[hii], mzpath+'sdss_mz_hii_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec[hii], mzpath+'sdss_mz_hii_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass[hii], mzpath+'sdss_mz_hii_mass.fits', clobber=clobber

       im_mwrfits, mz_parent[agn], mzpath+'sdss_mz_agn_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec[agn], mzpath+'sdss_mz_agn_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass[agn], mzpath+'sdss_mz_agn_mass.fits', clobber=clobber
    endif else begin            ; AGES
       splog, file=mzpath+'build_mz_emline.ages.log'
       splog, '#########################'
       splog, 'Building AGES emission-line sample'

; read the emission-line catalogs and subscript to the parent sample 
       ispec = read_ages(/ppxf)
       match, parent.ages_id, ispec.ages_id, m1, m2
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
       ispec = ispec[m2]

;      unfluxed = read_ages(/ppxf,/unfluxed)
;      unfluxed = unfluxed[m2]

; now select the emission-line sample; see MZPLOT_SAMPLE for details 
       splog, 'Flux cut = ', total(ispec.h_beta[0] gt ages_hbcut)
       
       mz = where($
         (ispec.h_beta[0] gt ages_hbcut) and $
;        (ispec.h_beta[0]/ispec.h_beta[1] gt ages_snrcut) and $ ; nominal cut
         (ispec.h_beta_ew[0] gt ages_ewcut) and $
         (ispec.oii_3727[1] ne -2.0) and $
         (ispec.oii_3727[0] gt 10^ages_oiihbcut*ispec.h_beta[0]),ngal)
;      mz = where($
;        (ispec.h_beta_ew[0] gt ewhbcut1) and $
;        (ispec.h_beta_ew[0]/ispec.h_beta_ew[1] gt ages_snrcut) and $
;        (ispec.oii_3727_ew[0]/ispec.oii_3727_ew[1] gt ages_snrcut) and $
;        (ispec.oiii_5007_ew[0]/ispec.oiii_5007_ew[1] gt ages_snrcut),ngal)
       splog, 'Emission-line sample = ', ngal

       mz_ispec = ispec[mz]
       mz_parent = parent[mz]
       mz_mass = mass[mz]

       splog, '[OIII] upper limits = ', total(mz_ispec.oiii_5007[1] eq -1), $
         total(mz_ispec.oiii_5007[1] eq -1)/float(ngal)
       
;      djs_plot, ispec.h_beta_ew[0], ispec.h_beta_ew[0]/ispec.h_beta_ew[1], $
;        psym=6, /xlog, /ylog, xr=[0.1,500], yr=[0.1,1000]
;      djs_oplot, mz_ispec.h_beta_ew[0], mz_ispec.h_beta_ew[0]/mz_ispec.h_beta_ew[1], $
;        psym=6, color='blue'
;      ww = where(mz_ispec.h_beta_ew[0] lt 1.0)
;      qaplot_ages_gandalf_specfit, mz_ispec[ww], psfile='junk.ps'

; classify into SF and AGN using various techniques
; candidate AGN at z>0.55
; * 201/183 - strong Ne III
; * 407/099 - strong Ne III

;      check = where((mz_ispec.nii_6584[1] eq -1.0) and (mz_ispec.h_alpha[1] eq -1.0))
;      qaplot_ages_gandalf_specfit, mz_ispec[check], psfile='junk.ps', /solar

; 1) broad lines 
       broad = (mz_ispec.isbroad eq 1) or (2.35*mz_ispec.h_beta_sigma[0] gt 800.0)
       splog, 'Broad ', total(broad), total(broad)/float(ngal)

; 2) [Ne V]
       nev = intarr(ngal)
       nev[ages_isnev(mz_ispec)] = 1
       splog, '[Ne V] ', total(nev), total(nev)/float(ngal)

; 3) BPT diagram, including upper limits; do not use the pure [NII]/Ha
; classification since we deal with upper limits below
       iclass = iclassification(mz_ispec,ratios=iratios1,$
         snrcut_class=0.0,silent=0,doplot=doplot,/noniihaclass)
;      class = doclass(mz_ispec,iratios=iratios,$
;        snrcut1=ages_snrcut,doplot=doplot)
;      bpt_agn = strmatch(iratios.final_class,'*AGN*')
       bpt_agn = (strtrim(iratios1.final_class,2) eq 'AGN') or $
         (strtrim(iratios1.final_class,2) eq 'SF/AGN')
       bpt_sf = strtrim(iratios1.final_class,2) eq 'SF'

       iratios = mzbpt_upper_limits(iratios1,mz_ispec,doplot=doplot)
       bpt_sf_limit = strmatch(iratios.final_class,'BPT_SF_*LIMIT*')
       bpt_agn_limit = strmatch(iratios.final_class,'BPT_AGN_*LIMIT*')

       splog, 'BPT-SF ', total(bpt_sf), total(bpt_sf)/float(ngal)
       splog, 'BPT-AGN ', total(bpt_agn), total(bpt_agn)/float(ngal)
       splog, 'BPT-SF-LIMITS ', total(bpt_sf_limit), total(bpt_sf_limit)/float(ngal)
       splog, 'BPT-AGN-LIMITS ', total(bpt_agn_limit), total(bpt_agn_limit)/float(ngal)

; note: a small fraction of objects (2080/198024=1.05%) cannot be
; classified using the above techniques either because their position
; in the BPT diagram was ambiguous or because H-alpha was not detected
; (1289/2080=62%)
       nobptclass = (strtrim(iratios.final_class,2) eq 'UNKNOWN')
       noha = (mz_ispec[where(nobptclass)].h_alpha[1] ne -2.0) and $
         (iratios[where(nobptclass)].nii_ha lt -900) and $
         (iratios[where(nobptclass)].nii_ha_limit lt -900)

       splog, 'No BPT classification = ', total(nobptclass), total(nobptclass)/float(ngal)
       splog, '  Of these, no H-alpha = ', total(noha), total(noha)/total(nobptclass)
       
; 4) Yan diagram; all the [OIII]/Hb upper limits are consistent with
; being star-forming
       ub = mz_parent.k_ubvrijhk_absmag_00[0]-mz_parent.k_ubvrijhk_absmag_00[1]
       yanagn = (iratios.oiii_hb ge ((1.4-1.2*ub)>(-0.1))) or $
         (iratios.oiii_hb_limit ge ((1.4-1.2*ub)>(-0.1)))
       yanandnobpt = (yanagn eq 1) and (nobptclass)
       yanandoiiilimit = (yanagn eq 1) and (iratios.oiii_hb_limit gt -900)

       splog, 'Yan AGN', total(yanagn), total(yanagn)/float(ngal)
       splog, '  Of these, no BPT class ', total(yanandnobpt), total(yanandnobpt)/float(ngal)
       splog, '  Of these, [OIII]/Hb limit ', total(yanandoiiilimit), total(yanandoiiilimit)/float(ngal)

; debugging plot
       if keyword_set(doplot) then begin
          det = where(iratios.oiii_hb gt -900.0,comp=lim)
          djs_plot, ub[det], iratios[det].oiii_hb, psym=6, xr=[0,2], sym=0.2
          djs_oplot, ub[lim], iratios[lim].oiii_hb_limit, psym=6, sym=0.2, color='yellow'
          ubaxis = range(0.0,2.0,500)
          djs_oplot, ubaxis, (1.4-1.2*ubaxis)>(-0.1), line=0, thick=7
          wdet = where(yanagn[det] eq 1)
          wlim = where(yanagn[lim] eq 1)
          djs_oplot, ub[det[wdet]], iratios[det[wdet]].oiii_hb, psym=6, sym=0.2, color='cyan'
;         djs_oplot, ub[wlim], iratios[wlim].oiii_hb_limit, psym=6, sym=0.2, color='yellow'
       endif
       
; 5) X-ray AGN; a handful of the sources with L(X)<10^41 and z<0.15
; are likely star-forming, but they represent such a small number of
; objects that we do not bias our sample by excluding them
       xray = mz_parent.x_match eq 1
       splog, 'Xray ', total(xray), total(xray)/float(ngal)

; 6) IRAC AGN
       irac = ages_irac_agn(mz_parent,debug=debug)
       splog, 'IRAC ', total(irac), total(irac)/float(ngal)

; 7) radio AGN; convert from the radio flux at 1.4GHz in mJy to Power
; in W/Hz (see Hickox+09); assume a power-law intrinsic spectrum to
; compute the K-correction (see pg 166 in Peterson's book) and adopt
; the AGN cut from Kauffmann+08 (see also Ryan's paper)
       radio = intarr(ngal)
       isradio = where(mz_parent.wsrt_match,nisradio)
       if (nisradio ne 0L) then begin
          alpha = 0.5
          kcorr = (1+mz_parent[isradio].z)^(alpha-1) ; eq 10.29 in Peterson
          pradio = mz_parent[isradio].wsrt_flux*kcorr*1E-3*1D-26*$ ; [W/Hz]
            4.0*!dpi*dluminosity(mz_parent[isradio].z,/meter)^2
;         plot, mz_parent[isradio].z, pradio, psym=6, /ylog
;         oplot, !x.crange, 10^23.7*[1,1]
          radio[isradio] = alog10(pradio) gt 23.8
       endif
       splog, 'Radio ', total(radio), total(radio)/float(ngal)

; pack it all into a structure
       class = replicate({agn: 0, broad_agn: 0, nev_agn: 0, $
         bpt_agn: -1, yan_agn: 0, xray_agn: 0, irac_agn: 0, radio_agn: 0},ngal)
       class = struct_addtags(class,im_struct_trimtags(iratios,$
         select=['bpt_d','bpt_phi','nii_ha','nii_ha_err','oiii_hb',$
         'oiii_hb_err','nii_ha_limit','oiii_hb_limit'],$
         newtags=['bpt_d','bpt_phi','bpt_nii_ha','bpt_nii_ha_err','bpt_oiii_hb',$
         'bpt_oiii_hb_err','bpt_nii_ha_limit','bpt_oiii_hb_limit']))

       class[where(bpt_agn or bpt_agn_limit)].bpt_agn = 1 ; 1=agn, 0=sf, -1=unclassified
       class[where(bpt_sf or bpt_sf_limit)].bpt_agn = 0
       
       class.broad_agn = broad
       class.nev_agn = nev
       class.yan_agn = yanagn
       class.xray_agn = xray
       class.irac_agn = irac
       class.radio_agn = radio
       
; final classifications: final class is AGN if it's been classified as
; an AGN using any of the other criteria above, unless it's been
; classified as star-forming using the BPT diagram
       class.agn = (class.bpt_agn eq 1) or (class.broad_agn eq 1) or $
         (class.nev_agn eq 1) or (class.yan_agn eq 1) or $
         (class.xray_agn eq 1) or (class.irac_agn eq 1) or $
         (class.radio_agn eq 1)
       reset = where(class.agn eq 1 and class.bpt_agn eq 0,nreset)
       if (nreset ne 0L) then class[reset].agn = 0
       
; some statistics
;      nclass = total(class.bpt_agn eq -1)
;      yclass = total(class.bpt_agn ne -1)
;      yagn = total((class.bpt_agn ne -1) and (class.bpt_agn eq 1))
;      ysf = total((class.bpt_agn ne -1) and (class.bpt_agn eq 0))
       splog, 'Final statistics:'
       splog, 'AGN: ', total(class.agn eq 1), ngal, total(class.agn eq 1)/float(ngal)
       splog, 'SF : ', total(class.agn eq 0), ngal, total(class.agn eq 0)/float(ngal)
;      im_plothist, mz_ispec[where(class.bpt_agn ne -1)].z, bin=0.02, xr=[0,0.8]
;      im_plothist, mz_ispec[where(class.bpt_agn eq -1)].z, bin=0.02, /over, line=5, /fill

; store the final classifications and then write out 
       mz_ispec = struct_addtags(temporary(mz_ispec),class)
       hii = where(mz_ispec.agn eq 0,nhii,comp=agn,ncomp=nagn)

       im_mwrfits, mz_parent, mzpath+'ages_mz_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec, mzpath+'ages_mz_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass, mzpath+'ages_mz_mass.fits', clobber=clobber

       im_mwrfits, mz_parent[hii], mzpath+'ages_mz_hii_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec[hii], mzpath+'ages_mz_hii_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass[hii], mzpath+'ages_mz_hii_mass.fits', clobber=clobber

       im_mwrfits, mz_parent[agn], mzpath+'ages_mz_agn_ancillary.fits', clobber=clobber
       im_mwrfits, mz_ispec[agn], mzpath+'ages_mz_agn_ispec.fits', clobber=clobber
       im_mwrfits, mz_mass[agn], mzpath+'ages_mz_agn_mass.fits', clobber=clobber
    endelse
    splog, /close
    
return
end


;; Lama diagram    
;      xx = alog10(mz_ispec.oii_3727_ew[0]/mz_ispec.h_beta_ew[0])
;      yy = alog10(mz_ispec.oiii_5007_ew[0]/mz_ispec.h_beta_ew[0])
;      oiiihblamaagn = 0.14/(xx-1.45)+0.83
;      lamaagn = where(yy gt oiiihblamaagn,nlamaagn)
;      if (nlamaagn ne 0) then class[lamaagn].class = $
;        class[lamaagn].class+'/LAMA_AGN'

;       broad = where((mz_ispec.isbroad eq 1),nbroad)
;       if (nbroad ne 0L) then class[broad].class = $
;         class[broad].class+'/BROAD_AGN'
;
;       ntot = total(mz_ispec.isbroad eq 0)
;       isagn = (strmatch(iratios.final_class,'*AGN*') eq 1) and (mz_ispec.isbroad eq 0)
;       issf = (strtrim(class.class,2) eq 'SF') and (mz_ispec.isbroad eq 0)
;       splog, total(isagn), total(isagn)/ntot, total(issf), total(issf)/ntot
;
;       niifail = (iratios.nii_ha lt -900) and (mz_ispec.isbroad eq 0)
;       niifail_loz = (iratios.nii_ha lt -900) and (mz_ispec.z le 0.4) and $
;         (mz_ispec.isbroad eq 0)
;       niifail_hiz = (iratios.nii_ha lt -900) and (mz_ispec.z gt 0.4) and $
;         (mz_ispec.isbroad eq 0)
;       splog, total(niifail), total(niifail)/float(ngal), $
;         total(niifail_hiz), total(niifail_hiz)/total(niifail)
;
;       im_plothist, mz_parent[where(niifail)].k_ubvrijhk_absmag_00[1], bin=0.1
;       im_plothist, mz_parent[where(niifail_loz)].k_ubvrijhk_absmag_00[1], bin=0.1, /over, line=5
;       im_plothist, mz_parent[where(niifail_hiz)].k_ubvrijhk_absmag_00[1], bin=0.1, /over, line=3
       
;; [Ne V]
;       nevagn = where((mz_ispec.nev_3426[0]/mz_ispec.nev_3426[1] gt snrcut1),nnevagn)
;       ww = cmset_op(nevagn,'and',where(strmatch(class.class,'*AGN*') eq 0))
;       qaplot_ages_gandalf_specfit, mz_ispec[ww], psfile='junk.ps'
;       if (nnevagn ne 0) then class[nevagn].class = $
;         class[nevagn].class+'/NEV_AGN'

