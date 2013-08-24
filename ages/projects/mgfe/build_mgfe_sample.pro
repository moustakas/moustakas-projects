pro build_mgfe_sample
; jm13apr14siena

    mgfepath = getenv('IM_PROJECTS_DIR')+'/mgfe/'
    ppxf = read_ages(/ppxf)
    phot = read_ages(/phot)
    kcorr = read_ages(/kcorr)

    phot = phot[ppxf.ages_id]
    kcorr = kcorr[ppxf.ages_id]

    indx = lindgen(n_elements(phot))
;   indx = where(ppxf.vdisp ne 165.0 and ppxf.z gt 0.05 and ppxf.z lt 0.75)
;   indx = where(ppxf.oii_3727_ew[0] lt 1.0 and ppxf.d4000_narrow[0] gt 1.5 and $
;     ppxf.vdisp ne 165.0 and ppxf.z gt 0.05 and ppxf.z lt 0.7)
;   qaplot_ages_gandalf_specfit, ppxf[indx[0:10]], ss, psfile='~/junk.ps'

    specfit = read_ages_gandalf_specfit(ppxf[indx])

; write out
    im_mwrfits, specfit, mgfepath+'ages_mgfe_specfit.fits', /clobber
    im_mwrfits, ppxf[indx], mgfepath+'ages_mgfe_ppxf.fits', /clobber
    im_mwrfits, kcorr[indx], mgfepath+'ages_mgfe_kcorr.fits', /clobber
    
stop    
    
return
end
    
