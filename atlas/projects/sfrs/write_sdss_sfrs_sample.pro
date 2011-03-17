pro write_sdss_sfrs_sample, sdust, snodust, ancillary, makeplots=makeplots, write=write
; jm05aug24uofa - write an SDSS sample for the SFRs paper

    outpath = atlas_path(/projects)+'sfrs/'

    sigmacut = 3.0
    
    if (n_elements(sdust) eq 0L) then sdust = read_sdss(sdssnodust=snodust,ancillary=ancillary);,/dr2)

    cut = where((ancillary.z_obj gt 0.033) and (ancillary.z_obj lt 0.25) and $
      (ancillary.infiber_sdss_r gt 0.1) and (ancillary.infiber_sdss_r le 1.0),ncut)

    splog, string(ncut,format='(I0)')+'/'+string(n_elements(sdust),format='(I0)')+$
      ' galaxies made the fiber fraction and redshift cuts.'

    keep = where((sdust.h_alpha[0]/sdust.h_alpha[1] gt sigmacut) and $
      (sdust.h_beta[0]/sdust.h_beta[1] gt sigmacut) and $
      (ancillary.z_obj gt 0.033) and (ancillary.z_obj lt 0.25) and $
      (ancillary.infiber_sdss_r gt 0.1) and (ancillary.infiber_sdss_r le 1.0))

;   hakeep = where(sdust.h_alpha[0]/sdust.h_alpha[1] gt sigmacut,nha)
;   hbkeep = where(sdust.h_beta[0]/sdust.h_beta[1] gt sigmacut,nhb)
;   keep = cmset_op(cmset_op(hakeep,'and',hbkeep),'and',cut)

    toss = where((sdust.h_alpha[0]/sdust.h_alpha[1] le sigmacut) and $
      (sdust.h_beta[0]/sdust.h_beta[1] le sigmacut) and $
      (ancillary.z_obj gt 0.033) and (ancillary.z_obj lt 0.25) and $
      (ancillary.infiber_sdss_r gt 0.1) and (ancillary.infiber_sdss_r le 1.0))

    ngalaxy = n_elements(keep)

    splog, string(ngalaxy,format='(I0)')+'/'+string(n_elements(sdust),format='(I0)')+$
      ' galaxies made the S/N cuts.'

    hii = where(strtrim(sdust[keep].bpt_class,2) eq 'HII',nhii)
    agn = where(strtrim(sdust[keep].bpt_class,2) eq 'AGN',nagn)

    splog, 'Found '+string(nhii,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' star-forming galaxies.'
    splog, 'Found '+string(nagn,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' AGN.'
    
    if keyword_set(write) then begin

; --------------------       

       splog, 'Writing '+outpath+'sdss_sfrs_hii_speclinefit.fits.gz'
       mwrfits, sdust[keep[hii]], outpath+'sdss_sfrs_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'sdss_sfrs_hii_speclinefit_nodust.fits.gz'
       mwrfits, snodust[keep[hii]], outpath+'sdss_sfrs_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_hii_speclinefit_nodust.fits'], /sh

       splog, 'Writing '+outpath+'sdss_sfrs_hii_ancillary.fits.gz'
       mwrfits, ancillary[keep[hii]], outpath+'sdss_sfrs_hii_ancillary.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_hii_ancillary.fits'], /sh

; --------------------       
       
       splog, 'Writing '+outpath+'sdss_sfrs_agn_speclinefit.fits.gz'
       mwrfits, sdust[keep[agn]], outpath+'sdss_sfrs_agn_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_agn_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'sdss_sfrs_agn_speclinefit_nodust.fits.gz'
       mwrfits, snodust[keep[agn]], outpath+'sdss_sfrs_agn_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_agn_speclinefit_nodust.fits'], /sh

       splog, 'Writing '+outpath+'sdss_sfrs_agn_ancillary.fits.gz'
       mwrfits, ancillary[keep[agn]], outpath+'sdss_sfrs_agn_ancillary.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_sfrs_agn_ancillary.fits'], /sh

    endif

    if keyword_set(makeplots) then begin

; z/zmax       
       
       bin = 0.02
       plothist, ancillary.zzmax, bin=bin, xr=[0,1], thick=2.0, $
         charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
         xtitle=textoidl('z/z_{max}'), ytitle='Number', line=2, /nan
       plothist, ancillary[cut].zzmax, bin=bin, /overplot, /nan
       plothist, ancillary[keep].zzmax, bin=bin, /overplot, $
         /fill, /fline, color=fsc_color('light grey',100), fspacing=0.05, $
         orientation=45, /nan
       plothist, ancillary[toss].zzmax, bin=bin, /overplot, $
         /fill, color=fsc_color('black',90), fspacing=0.05, /nan

       cc = get_kbrd(1)
       
; apparent r magnitude
       
       bin = 0.1
       plothist, ancillary.sdss_r, bin=bin, xr=[13.5,18], thick=2.0, $
         charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
         xtitle=textoidl('r'), ytitle='Number', line=2, /nan
       plothist, ancillary[cut].sdss_r, bin=bin, /overplot, /nan
       plothist, ancillary[keep].sdss_r, bin=bin, /overplot, $
         /fill, /fline, color=fsc_color('light grey',100), fspacing=0.05, $
         orientation=45, /nan
       plothist, ancillary[toss].sdss_r, bin=bin, /overplot, $
         /fill, color=fsc_color('black',90), fspacing=0.05, /nan
       
       cc = get_kbrd(1)

; D(4000)
       
       bin = 0.05
       plothist, sdust.model_d4000_narrow[0], bin=bin, xr=[0.8,2.2], thick=2.0, $
         charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
         xtitle=textoidl('D(4000)'), ytitle='Number', line=2, /nan
       plothist, sdust[cut].model_d4000_narrow[0], bin=bin, /overplot, /nan
       plothist, sdust[keep].model_d4000_narrow[0], bin=bin, /overplot, $
         /fill, /fline, color=fsc_color('light grey',100), fspacing=0.05, $
         orientation=45, /nan
       plothist, sdust[toss].model_d4000_narrow[0], bin=bin, /overplot, $
         /fill, color=fsc_color('black',90), fspacing=0.05, /nan
       
       cc = get_kbrd(1)

    endif
    
return
end
    
