pro ages_refit_postbursts, refit=refit
; jm09feb19nyu - the emission-line measurements of Michael's
; post-starbursts are frequently crap; so refit them here

    base_specfitpath = ages_path(/specfit)
    base_spec1dpath = ages_path(/spec1d)
    version = ages_version(/ispec_specfit)
    spec1dpath = base_spec1dpath+'fluxed/after_skysubpca/'
    
    path = ages_path(/projects)+'postburst/brown/'
    sum1 = rsex(path+'summary.txt')

    aa1 = read_ages(/ancillary)
    ii1 = read_ages(/ispec)
    spherematch, aa1.ra, aa1.dec, 15.0*im_hms2dec(sum1.ra), $
      im_hms2dec(sum1.dec), 1.0/3600.0, m1, m2
    aa = aa1[m1]
    ii = ii1[m1]
    sum = sum1[m2]

    continuum_maxwave = 8500.0  ; exclude wavelengths contaminated by the red leak
    specres = 4.5               ; 5.0        ; FWHM instrumental resolution [A]
    vdisp = 50.0                ; assume constant stellar velocity dispersion [km/s]

    vlinemaxtweak = 500.0
    sigmax = 300.0

    indexfile = base_specfitpath+'indexlist_'+version+'.dat'
    linefile = base_specfitpath+'elinelist_'+version+'.dat'
    Zbc03 = 'Z02'
;   Zbc03 = ['Z004','Z02','Z05']    
    templatefile = base_specfitpath+'BC03_'+Zbc03+'_salpeter_ages_'+$
      ages_version(/templates)+'_templates.fits'

    allpass = aa.pass
    upass = allpass[uniq(allpass,sort(allpass))]
    npass = n_elements(upass)

    if keyword_set(refit) then begin
    
       t0 = systime(1)
;      for ii = 13, 13 do begin
       for ii = 0, npass-1 do begin

          case upass[ii] of
             105: continuum_minwave = 4600.0
             111: continuum_minwave = 4600.0
             112: continuum_minwave = 4500.0
             113: continuum_minwave = 4900.0
             115: continuum_minwave = 4700.0
             201: continuum_minwave = 4800.0
             203: continuum_minwave = 5300.0
             205: continuum_minwave = 4800.0
             206: continuum_minwave = 4300.0
             207: continuum_minwave = 4300.0
             208: continuum_minwave = 5400.0
             210: continuum_minwave = 4800.0
             213: continuum_minwave = 5100.0
             301: continuum_minwave = 5200.0
             303: continuum_minwave = 4600.0
             304: continuum_minwave = 5200.0
             307: continuum_minwave = 5200.0
             308: continuum_minwave = 5200.0
             309: continuum_minwave = 5200.0
             313: continuum_minwave = 4200.0
             314: continuum_minwave = 4500.0
             422: continuum_minwave = 4400.0
             else: delvarx, continuum_minwave
          endcase             
          
          ww = where(upass[ii] eq allpass,nww)
          suffix = 'refit_'+string(upass[ii],format='(I3.3)')
          thispass = string(aa[ww[0]].pass,format='(I3.3)')
          data = mrdfits(spec1dpath+'ages_'+thispass+'.fits.gz',1,/silent)

          index = lonarr(nww)
          for jj = 0L, nww-1L do index[jj] = where(aa[ww[jj]].aper eq data.aper)

          galaxy = data.galaxy[index]
          zobj = data.z[index]
          wave = data.wave
          flux = data.flux
          invvar = (1.0D/(data.ferr+(data.ferr eq 0.0))^2.0)*(data.ferr ne 0.0)

          case upass[ii] of
             302: zsnrcross = -1.0
             403: zsnrcross = -1.0
             else: zsnrcross = 2.0
          endcase
          
          specdata = ispeclinefit(wave[*,index],flux[*,index],invvar[*,index],$
            zobj=zobj,specres=specres,vdisp=vdisp,sigmax=sigmax,$
            zsnrcross=zsnrcross,galaxy=galaxy,outpath=specfitpath,$
            suffix=suffix,templatefile=templatefile,linefile=linefile,indexfile=indexfile,$
            specfit=specfit,continuum_minwave=continuum_minwave,$
            continuum_maxwave=continuum_maxwave,/nologfile,vmaxtweak=vmaxtweak,$
            vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,doplot=doplot,$
            specdatafile=specdatafile)

       endfor
       splog, 'Total time = ', (systime(1)-t0)/60.0

    endif

; merge and parse    
    
    allfits = file_search(path+'?????_refit_???_specdata.fits.gz',count=nn)
    for kk = 0L, nn-1L do begin
       fits1 = mrdfits(allfits[kk],1,/silent)
       if (kk eq 0L) then fits = fits1 else fits = [fits,fits1]
    endfor

    occ = iclassification(ii,doplot=0,snrcut_class=2.0,ratios=orr)
    cc = iclassification(fits,doplot=0,snrcut_class=2.0,ratios=rr)
;   struct_print, rr
    final = struct_trimtags(rr,select=['galaxy','nii_ha*','oiii_hb*',$
      'final_class'])
    
; in the following objects ICLASS is wrong - there are no lines 
    ww = speclinefit_locate(final,'ages_'+['112/001','204/105'])
    final[ww].nii_ha = -999.0 & final[ww].nii_ha_err = -999.0
    final[ww].oiii_hb = -999.0 & final[ww].oiii_hb_err = -999.0
    final[ww].final_class = 'UNKNOWN'

; in the following compute a limit on [NII]/Ha
    ww = speclinefit_locate(final,'ages_'+['111/067','115/232','207/160'])
    final[ww].nii_ha = alog10(fits[ww].nii_6584_limit/fits[ww].h_alpha[0])
    final[ww].nii_ha_err = -1.0
    final[ww[0]].final_class = 'SF'
    final[ww[1]].final_class = 'SF/AGN' ; 'SF_or_SF/AGN'
    final[ww[2]].final_class = 'SF'

    openw, lun, path+'postburst_class.dat', /get_lun
    printf, lun, '## Optical classifications for the 24 AGES post-starbursts in Brown et al. (2009)'
    printf, lun, '## J. Moustakas, NYU'
    printf, lun, '## '+im_today()
    printf, lun, '## '
    printf, lun, '## SF=below the Kauffmann+03 curve'
    printf, lun, '## AGN=above the Kewley+01 curve'
    printf, lun, '## SF/AGN=between the two curves'
    printf, lun, '## UNKNOWN=no classification possible (2 or more lines not detected)'
    printf, lun, '## '
    printf, lun, '## Notes: All detected lines have S/N>2; undetected lines have '
    printf, lun, '##   line-ratios and errors equal to -999.  In three objects [NII] '
    printf, lun, '##   was not detected; here, the quoted [NII]/Ha ratio is a 1-sigma '  
    printf, lun, '##   upper limit and the error on [NII]/Ha is set to -1.0.  One of ' 
    printf, lun, '##   these objects, ages_115/232, formally lies in the SF/AGN region '
    printf, lun, '##   of the diagram, but it could also be in the SF region.'
    printf, lun, '## '
    printf, lun, '# 1 GALAXY'
    printf, lun, '# 2 NII_HA [log]'
    printf, lun, '# 3 NII_HA_ERR [log]'
    printf, lun, '# 4 OIII_HB [log]'
    printf, lun, '# 5 OIII_HB_ERR [log]'
    printf, lun, '# 6 CLASS'
    struct_print, final, lun=lun, /no_head
    free_lun, lun
    
stop    
    
return
end
    
    
