pro compare_ispec_platefit, sdssancillary1, sdssdust1, spec1d=spec1d, specfit=specfit
; jm07sep11nyu - compare the emission-line fluxes and EWs measured
;                from spectra that have been fitted using both ispec
;                and C. Tremonti's platefit

; some datapaths and fitting parameters    
    
    datapath = ages_path()+'compare_ispec_platefit/'
    spec1dpath = datapath+'spec1d/'

    eigendir = ages_path(/templates)
    eigenfile = 'BC03_Z02_salpeter_templates.fits'
    
    medsmooth_window = 250L     ; larger than the default values!
    smooth_window = 100L

; read the SDSS database and pick some spectra to fit    
    
    if (n_elements(sdssancillary1) eq 0L) then sdssancillary1 = read_sdss(specdata=sdssdust1)

    allplates = sdssancillary1[uniq(sdssancillary1.plateid,sort(sdssancillary1.plateid))].plateid
    plate = allplates[400]
    these = where(sdssancillary1.plateid eq plate)
    ngalaxy = n_elements(these)

    sdssancillary = sdssancillary1[these]
    sdssdust = sdssdust1[these]

    mjd = (sdssancillary.mjd)[0] ; should be the same for all fibers

; read the spectra themselves and pack everything into a structure     
    
    splog, 'Reading '+string(ngalaxy,format='(I0)')+' SDSS spectra.'
    readspec, plate, sdssancillary.fiberid, mjd=mjd, flux=flux, flerr=flerr, $
      invvar=invvar, disp=disp, wave=wave, /align

    specpass = {sdss_id: sdssancillary.sdss_id, ra: sdssancillary.ra, $
      dec: sdssancillary.dec, plate: sdssancillary.plateid, $
      fiber: sdssancillary.plateid, mjd: sdssancillary.mjd, $
      z: sdssancillary.z, vdisp: sdssancillary.v_disp, $
      wave: wave, flux: flux*1D-17, ferr: flerr*1D-17}

; now do the fitting

    if keyword_set(specfit) then begin

       specdata = sdss_ispeclinefit(specpass,specres=2.0,snrcut=0.0,dustmodel=0,zsnrcross=2.0,$
         datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,sigmax=800.0,$
         /charlot,/zcrosscor,postscript=1,/write,vmaxshift=400.0,nologfile=1,$
         starvdisp=starvdisp,doplot=0,debug=0,eigenfile=eigenfile,eigendir=eigendir,$
         medsmooth_window=medsmooth_window,smooth_window=smooth_window)

    endif

    stop    

    w = where(specdata.h_beta_ew[0]/specdata.h_beta_ew[1] gt 3.0)
    plot, specdata[w].h_beta_ew[0], 100*(specdata[w].h_beta_ew[0]-sdssdust[w].h_beta_ew[0])/specdata[w].h_beta_ew[0], ps=4, xr=[-3,25], yr=[-100,100]
    jj = im_stats(100*(specdata[w].h_beta_ew[0]-sdssdust[w].h_beta_ew[0])/specdata[w].h_beta_ew[0],/ver,sigrej=3.0)

    w = where(specdata.oii_3727_ew[0]/specdata.oii_3727_ew[1] gt 3.0)
    plot, specdata[w].oii_3727_ew[0], 100*(specdata[w].oii_3727_ew[0]-sdssdust[w].oii_3727_ew[0])/specdata[w].oii_3727_ew[0], ps=4, xr=[-3,25], yr=[-100,100]
    jj = im_stats(100*(specdata[w].oii_3727_ew[0]-sdssdust[w].oii_3727_ew[0])/specdata[w].oii_3727_ew[0],/ver,sigrej=3.0)

    w = where(specdata.oiii_5007_ew[0]/specdata.oiii_5007_ew[1] gt 3.0)
    plot, specdata[w].oiii_5007_ew[0], 100*(specdata[w].oiii_5007_ew[0]-sdssdust[w].oiii_5007_ew[0])/specdata[w].oiii_5007_ew[0], ps=4, xr=[-3,25], yr=[-100,100]
    jj = im_stats(100*(specdata[w].oiii_5007_ew[0]-sdssdust[w].oiii_5007_ew[0])/specdata[w].oiii_5007_ew[0],/ver,sigrej=3.0)



return
end
    
