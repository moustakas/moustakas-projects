pro write_00jansen, data, writeancillary=writeancillary
; parse the NFGS integrated spectrophotometric measurements 

; jm03sep10uofa - originally written
; jm04dec10uofa - major re-write    

; NOTE: E(B-V) [BH] = A_B/1.324/3.1 ; Jansen, private communication    
    
; the following objects are not part of the statistical sample: 8, 19
    
    root = '00jansen'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    filt = im_filterspecs(filternames=['bessell_U.par',$
      'bessell_B.par','bessell_V.par','bessell_R.par'])
    
; ---------------------------------------------------------------------------    
; preliminaries    
; ---------------------------------------------------------------------------    

    if keyword_set(writeancillary) then begin

       parse_ned_byname, '00jansen.ned', ned, inputnedfile='galaxy_list.txt', $
         outfile='00jansen_ned.fits', nedpath=path, outpath=path
       parse_ned_photometry, '00jansen.photometry.ned', photo, inputnedfile='galaxy_list.txt', $
         outfile='00jansen_ned_photo.fits', nedpath=path, outpath=path
       
       basicname = '00jansen_ned.fits.gz'
       photoname = '00jansen_ned_photo.fits.gz'
       distname = '00jansen_distances.fits.gz'
       diamname = '00jansen_diameters.fits.gz'

; don't RA-sort because NGC5407 (147) and NGC5425 (148) switch
; positions from Rolf's ordering; maintain the native NFGS ordering 
       
       write_ancillary_data, ancillary, datapath=path, outpath=outpath, $
         basicname=basicname, photoname=photoname, distname=distname, $
         diamname=diamname, /norasort 
       outname = '00jansen_ancillary_data.fits'
       splog, 'Writing '+path+outname+'.'
       mwrfits, ancillary, path+outname, /create
       spawn, ['gzip -f '+path+outname], /sh

    endif

; ---------------------------------------------------------------------------    

; read the data sample structure
    
    sample = rsex(path+'nfgs_sample_properties.dat')
    samplephoto = rsex(path+'nfgs_synth_photometry.dat')
    samplephoto2 = rsex(path+'photometry_table3.dat')

    sample = im_struct_trimtags(sample,select=tag_names(sample),newtags='NFGS_'+tag_names(sample))
    samplephoto = im_struct_trimtags(samplephoto,select=tag_names(samplephoto),newtags='NFGS_'+tag_names(samplephoto))
    ngalaxy = n_elements(sample)
    
    ancillary = mrdfits(path+root+'_ancillary_data.fits.gz',1,/silent)
    ancillary = struct_addtags(replicate({nfgs_galaxy: '', nfgs_ugc: ''},ngalaxy),ancillary)

    ancillary.nfgs_galaxy = sample.nfgs_name
    ancillary.nfgs_ugc = sample.nfgs_ugc
    sample = struct_trimtags(sample,except=['NFGS_NAME','NFGS_UGC'])

    data = init_cat_linefit(ngalaxy=ngalaxy)
    data = struct_trimtags(data,except=['GALAXY','ALT_GALAXY','RA','DEC','Z_OBJ','Z_OBJ_ERR',$
      'DISTANCE','DISTANCE_ERR','M_B','M_B_ERR','B_LUM','B_LUM_ERR'])

    data = struct_addtags(ancillary,struct_addtags(sample,data))
;   niceprint, data.nfgs_id, data.nfgs_galaxy, data.galaxy, data.nfgs_ugc, data.ned_galaxy
;   niceprint, data.nfgs_id, data.galaxy, ancillary.galaxy, ancillary.ned_galaxy

    linename = data[0].linename
    nline = n_elements(linename)

; ---------------------------------------------------------------------------    
; photometric properties from Jansen et al.; assume that BTOT and
; BVSPEC from Jansen et al. have been corrected for foreground
; Galactic extinction using Burstein & Heiles; undo this correction,
; and then apply the Schlegel et al. extinction correction 
; ---------------------------------------------------------------------------    

    AB_bh84 = sample.nfgs_a_b
    AV_bh84 = 0.755*AB_bh84 ; using k_lambda()
    AU_bh84 = 1.156*AB_bh84 ; from Jansen 
    AR_bh84 = 0.588*AB_bh84 ; from Jansen 

    AB_sfd = ancillary.ebv_mw*k_lambda(filt[1].weff,/odonnell)
    AV_sfd = ancillary.ebv_mw*k_lambda(filt[2].weff,/odonnell)
    AU_sfd = ancillary.ebv_mw*k_lambda(filt[0].weff,/odonnell)
    AR_sfd = ancillary.ebv_mw*k_lambda(filt[3].weff,/odonnell)
    
    good = where((samplephoto.nfgs_btot gt -900.0))
    data[good].b_err = samplephoto[good].nfgs_ebtot
    data[good].b = (samplephoto[good].nfgs_btot + AB_bh84[good]) - AB_sfd[good]

    good = where((samplephoto.nfgs_bvspec gt -900.0) and (samplephoto.nfgs_btot gt -900))
    data[good].v_err = sqrt(samplephoto[good].nfgs_ebtot^2 + 0.065^2)
    data[good].v = ((samplephoto[good].nfgs_btot - samplephoto[good].nfgs_bvspec) + AV_bh84[good]) - AV_sfd[good]

;   data[good].v = samplephoto[good].nfgs_btot - samplephoto[good].nfgs_bvspec + $
;     (AB_bh84[good]-AV_bh84[good]) - (AB_sfd[good]-AV_sfd[good])
    
    good = where((samplephoto2.ube gt -900.0) and (samplephoto.nfgs_btot gt -900))
    data[good].u_err = sqrt(samplephoto2[good].eube^2 + samplephoto2[good].ebtot^2)
;   data[good].u = samplephoto2[good].btot + samplephoto2[good].ube + $
;     (AU_bh84[good]-AB_bh84[good]) - (AU_sfd[good]-AB_sfd[good])
    data[good].u = ((samplephoto[good].nfgs_btot + samplephoto2[good].ube) + AU_bh84[good]) - AU_sfd[good]

    good = where((samplephoto2.bre gt -900.0) and (samplephoto.nfgs_btot gt -900))
    data[good].r_err = sqrt(samplephoto2[good].ebre^2 + samplephoto2[good].ebtot^2)
    data[good].r = ((samplephoto[good].nfgs_btot - samplephoto2[good].bre) + AR_bh84[good]) - AR_sfd[good]

;   niceprint, data.nfgs_galaxy, data.b, data.b_err, data.rc3_b, data.rc3_b_err
;   niceprint, data.u, data.u_err, data.b, data.b_err, data.v, data.v_err

    splog, 'Computing absolute magnitudes and luminosities.'

    band = ['U','B','V','R']
    tags = ['u','b','v','r']
    tags_err = tags+'_err'
    abstags = ['m_u','m_b','m_v','m_r']
    abstags_err = abstags+'_err'
    
    lumtags = ['u_lum','b_lum','v_lum','r_lum']
    lumtags_err = lumtags+'_err'

    for iband = 0L, n_elements(tags)-1L do begin

       true = tag_exist(data,tags[iband],index=tagsindx)
       true = tag_exist(data,tags_err[iband],index=tagsindx_err)
       true = tag_exist(data,abstags[iband],index=abstagsindx)
       true = tag_exist(data,abstags_err[iband],index=abstagsindx_err)
       true = tag_exist(data,lumtags[iband],index=lumtagsindx)
       true = tag_exist(data,lumtags_err[iband],index=lumtagsindx_err)

       good = where((data.distance gt -900.0) and (data.(tagsindx) gt -900.0),ngood)
       if (ngood ne 0L) then begin

          mags = im_absolute_magnitudes(band[iband],data[good].(tagsindx),$
            data[good].distance,mag_err=data[good].(tagsindx_err),$
            distance_err=data[good].distance_err)

          data[good].(abstagsindx)     = mags.absmag
          data[good].(abstagsindx_err) = mags.absmag_err
          
          data[good].(lumtagsindx)     = mags.lum
          data[good].(lumtagsindx_err) = mags.lum_err
          
       endif
       
    endfor

; ---------------------------------------------------------------------------
; Emission-Line Fluxes; append H-beta
; ---------------------------------------------------------------------------

    balmerfluxes = rsex(path+'nfgs_balmer_fluxes.dat')
    fluxes = rsex(path+'nfgs_fluxes.dat')
    fluxes = struct_addtags(fluxes,replicate({h_beta: 0.0},ngalaxy))
;   fluxes = struct_addtags(fluxes,im_struct_trimtags(balmerfluxes,$
;     select='H_BETA_ABS',newtags='H_BETA'))

    balmerews = rsex(path+'nfgs_balmer_EWs.dat')
    ews = rsex(path+'nfgs_EWs.dat')
;   ews = struct_addtags(ews,im_struct_trimtags(balmerews,$
;    select='H_BETA_EW_ABS',newtags='H_BETA'))

; figure out if the emission lines were normalized by H-beta or
; H-alpha
    
    norm = fltarr(ngalaxy)
    hanorm = where(fluxes.h_alpha eq 1.0,comp=hbnorm)
    norm[hanorm] = balmerfluxes[hanorm].h_alpha_raw
    norm[hbnorm] = balmerfluxes[hbnorm].h_beta_raw

; now fill the emission-line structure

    tags = tag_names(data)
    fluxtags = tag_names(fluxes)
    ewtags = tag_names(ews)

    for iline = 0L, nline-1L do begin

       match = where(strmatch(tags,linename[iline]) eq 1B,nmatch)

       wavematch = where(strmatch(tags,linename[iline]+'_WAVE') eq 1B,nwavematch)
       refmatch = where(strmatch(fluxtags,linename[iline]) eq 1B,nrefmatch)

       ewmatch = where(strmatch(tags,linename[iline]+'_EW') eq 1B,newmatch)
       refewmatch = where(strmatch(ewtags,linename[iline]) eq 1B,nrefewmatch)
          
       if (nrefmatch ne 0L) then begin
          
; H-BETA       

          if (strmatch(linename[iline],'*H_BETA*') eq 1B) then begin

             good = where((balmerfluxes.h_beta_abs gt 0.0) and (norm gt 0.0),ngood) ; flux, reddening-corrected             
             if (ngood ne 0L) then begin 

                wave = data[good].(wavematch)

                flux = balmerfluxes[good].h_beta_abs * 10^(0.4*k_lambda(wave,/odonnell)*data[good].ebv_mw)
                ferr = flux*0.15
                
                data[good].(match) = transpose([ [flux], [ferr] ])

; ###########################################################################                
                splog, 'WARNING: PROBLEM HERE!!!'
; ###########################################################################                

             endif
             
             good = where((balmerews.h_beta_ew_abs ne 0.0),ngood) ; EW             
             if (ngood ne 0L) then begin 

                ew = -balmerews[good].h_beta_ew_abs
                ewerr = ew*0.15
                
                data[good].(ewmatch) = transpose([ [ew], [ewerr] ])

             endif

             continue

          endif
          
; H-ALPHA

          if (strmatch(linename[iline],'*H_ALPHA*') eq 1B) then begin

             good = where((balmerfluxes.h_alpha_abs gt 0.0) and (norm gt 0.0),ngood) ; flux, reddening-corrected             
             if (ngood ne 0L) then begin 

                wave = data[good].(wavematch)

                flux = balmerfluxes[good].h_alpha_abs * 10^(0.4*k_lambda(wave,/odonnell)*data[good].ebv_mw)
                ferr = flux*0.15
                
                data[good].(match) = transpose([ [flux], [ferr] ])

             endif

             good = where((balmerews.h_alpha_ew_abs ne 0.0),ngood) ; EW
             if (ngood ne 0L) then begin
                
                ew = -balmerews[good].h_alpha_ew_abs
                ewerr = ew*0.15
                
                data[good].(ewmatch) = transpose([ [ew], [ewerr] ])

             endif

             continue
             
          endif
          
          good = where((fluxes.(refmatch) gt 0.0) and (norm gt 0.0),ngood) ; flux, reddening-corrected             
          if (ngood ne 0L) then begin

             wave = data[good].(wavematch)

             flux = norm[good]*fluxes[good].(refmatch) * 10^(0.4*k_lambda(wave,/odonnell)*data[good].ebv_mw)
             ferr = flux*0.15
             
             data[good].(match) = transpose([ [flux], [ferr] ])

          endif
          
          good = where((ews.(refewmatch) gt 0.0),ngood) ; EW
          if (ngood ne 0L) then begin
             
             ew = ews[good].(refewmatch)
             ewerr = ew*0.15
             
             data[good].(ewmatch) = transpose([ [ew], [ewerr] ])

          endif

       endif ; close NMATCH IF statement
          
    endfor

; ---------------------------------------------------------------------------
; Write out
; ---------------------------------------------------------------------------

; remove 2 objects that are not part of the statistical sample

;   indx = lindgen(ngalaxy)
;   remove, [[8,19]], indx
;   data = data[indx]

; de-redden the fluxes and then compute abundances

    idata = iunred_linedust(data,snrcut=3.0,/silent)
    adata = im_abundance(idata,snrcut=3.0)

    data = struct_addtags(data,adata)
    
    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; reddening-corrected fluxes    
    
    idata = struct_addtags(idata,adata)
    
    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, idata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh
    
return
end 
