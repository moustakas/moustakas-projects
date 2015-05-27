;+
; NAME:
;   BUILD_DESI_BGS_TEMPLATES
;
; PURPOSE:
;   Generate emission-line galaxy spectra for DESI using the DEEP2
;   spectroscopy and photometry.
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   minwave - minimum output wavelength (default 3000 A)
;   maxwave - approximate output wavelength (default 10,000 A) 
;   velpixsize - spectrum pixel size [default 20 km/s]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   I have verified that my BUILD_EMLINE script yields the same
;   emission-line spectrum as GANDALF, which was used to fit the DEEP2
;   spectra. 
; 
;   All line fluxes are defined with respect to [OII] 3729.
; 
; ToDo:
;   1. Add more nebular emission lines.
;   2. Allow the nebular emission-line ratios to vary.
;   3. Build stacks of spectra.
;   4. Add some jitter to the emission-line redshifts compared to the
;      continuum redshifts.  Could also vary the line-widths somewhat.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Jan 09, Siena
;   jm14mar11siena - major update
;
; Copyright (C) 2014, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function build_emline, emspectrum, logwave=logwave, lineflux=lineflux, $
  linesigma=linesigma, linewave=linewave, zshift=zshift
; internal support routine: build a single emission line
;   logwave - desired (rest-frame) wavelength vector [A]
;   lineflux - integrated (rest-frame) line-flux [erg/s/cm2]
;   linesigma - kinematic velocity width [km/s]
;   linewave - rest-frame emission-line wavelength [A]
;   zshift - radial velocity shift relative to the systemic redshift 

; allow the line to have a radial velocity shift with respect to the
; systemic redshift/velocity
    if n_elements(zshift) eq 0 then zshift = 0.0 

    light = 2.99792458D5
    sigma = linesigma/light/alog(10)       ; line-width [log-10 Angstrom]
    amplitude = lineflux/alog(10)/linewave ; line-amplitude [erg/s/cm2/A]
    linewave1 = alog10(linewave*(1+zshift))
    emspectrum += amplitude*exp(-0.5*(logwave-linewave1)^2/$ ; [erg/s/cm2/A, rest]
      sigma^2)/(sqrt(2.0*!pi)*sigma)

return, emspectrum
end

function continuum_convolve, in, velscale=velscale, vdisp=vdisp
; convolve to the AGES instrumental resolution
    fsps_inst_vdisp = 70.0 ; [km/s] (assume MILES models)
    thisvdisp = sqrt((vdisp^2.0-fsps_inst_vdisp^2)>0.0)
    if thisvdisp gt 20.0 then begin
       smoothing = thisvdisp/velscale ; [pixel]
       nkpix = long(4.0*ceil(smoothing))*2L+3
       klam = findgen(nkpix)-float(nkpix-1.0)/2.0
       kernel = exp(-0.5D*(klam/smoothing)^2)/sqrt(2.0*!dpi)/smoothing
       kernel = kernel/total(kernel)
       out = convol(in,kernel,/edge_truncate)
    endif else out = in
return, out
end

pro build_desi_bgs_templates, match_sdss=match_sdss, debug=debug, clobber=clobber
; jm15apr24siena - build a set of full-resolution templates for the
; BGS based on the AGES galaxy sample

;   echo "build_desi_bgs_templates" | /usr/bin/nohup idl > & ~/bgs-templates.log & 

    version = desi_bgs_templates_version()
    isedfit_version = desi_bgs_templates_version(/isedfit)

    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/desi/templates/'+$
      'bgs_templates/isedfit/'+isedfit_version+'/'
    templatepath = getenv('DESI_ROOT')+'/spectro/templates/bgs_templates/'+version+'/'
    sdssdr12dir = '/home/work/data/sdss/dr12/'
    isedfit_paramfile = isedfit_dir+'desi_ages_paramfile.par'

    if file_test(templatepath,/dir) eq 0 then file_mkdir, templatepath

    outfile_obs = templatepath+'bgs_templates_obs_'+version+'.fits'
    outfile_rest = templatepath+'bgs_templates_'+version+'.fits'
    outfile_continuum_rest = templatepath+'bgs_continuum_templates_'+version+'.fits'
    
; ---------------------------------------------------------------------------
; match against the SDSS
    if keyword_set(match_sdss) then begin
       sdssdir = getenv('IM_ARCHIVE_DIR')+'/projects/desi/templates/'+$
         'bgs_templates/sdss/'
       phot1 = mrdfits(sdssdir+'bootes-sdss-galaxies-dr12.fits',1)
       
       cat = mrdfits(outfile_obs,1)
       ngal = n_elements(cat)
       
       spherematch, phot1.ra, phot1.dec, cat.ra, cat.dec, 1D/3600.0, m1, m2
       srt = sort(m2)
       m1 = m1[srt] & m2 = m2[srt]

       phot = im_empty_structure(phot1,ncopies=ngal,empty_value=-999.0)
       phot[m2] = phot1[m1]

; try matching the missing objects against the stars
       phot1 = mrdfits(sdssdir+'bootes-sdss-stars-dr12.fits',1)

       miss = where(phot.run eq -999)
       spherematch, phot1.ra, phot1.dec, cat[miss].ra, cat[miss].dec, 1D/3600.0, m1, m2

       phot[miss[m2]] = im_struct_assign(phot1[m1],phot[miss[m2]])
       
; match against the spectroscopic catalog       
       spec1 = mrdfits(sdssdr12dir+'specObj-dr12.fits',1)
       spherematch, spec1.plug_ra, spec1.plug_dec, cat.ra, cat.dec, 1D/3600.0, m1, m2

       spec = im_empty_structure(spec1,ncopies=ngal,empty_value=-999.0)
       spec[m2] = spec1[m1]

       sdss_outfile = templatepath+'bgs_sdss_'+version+'.fits'
       mwrfits, phot, sdss_outfile, /create
       mwrfits, spec, sdss_outfile       

; write out some spectra
       jj = mrdfits('bgs_templates_sdss_v1.0.fits',2)
;      gg = (where(jj.fiberid gt -900))[1]
       gg = where(jj.plate eq 3864)
       niceprint, jj[gg].fiberid, jj[gg].plate, jj[gg].mjd
       readspec, jj[gg].plate, jj[gg].fiberid, mjd=jj[gg].mjd, flux=flux, wave=wave, zans=zans
       return
    endif

    if im_file_test(outfile_obs,clobber=clobber) then return
    if im_file_test(outfile_rest,clobber=clobber) then return
    if im_file_test(outfile_continuum_rest,clobber=clobber) then return
    
    params = read_isedfit_paramfile(isedfit_paramfile)
    weff = k_lambda_eff(filterlist=strtrim(params.filterlist,2))
    airtovac, weff, weff_vac

    light = im_light(/kms)
    light_ang = im_light(/Ang)
    maggies2flam = 10^(-0.4*48.6)*light_ang/weff^2
    
    synth_filters = (decam_filterlist())[[1,2,4]] ; =grz

    sfhgrid = 2 ; SFHGRID02 includes emission lines
    
    oiisnrcut = 3.0
    oiiewcut = 1.0

    rlimit = 23.5
    pc10 = 3.085678D19          ; fiducial distance [10 pc in cm]

; build the output *rest-wavelength* array; sample the rest-frame UV
; and optical wavelengths at high resolution, and the near-infrared at
; lower resolution (since it will only be used for filter
; convolutions) 
    
; simulation parameter defaults (rest-frame *vacuum* wavelengths!);
    velpixsize_hires = 20D ; [km/s]
    pixsize_hires = velpixsize_hires/light/alog(10) ; [pixel size in log-10 A]
    minwave_hires = alog10(1250D)
    maxwave_hires = alog10(4D4)
;   maxwave_hires = alog10(10000D)
    npix_hires = round((maxwave_hires-minwave_hires)/pixsize_hires+1L)

    restwave_hires = minwave_hires+dindgen(npix_hires)*pixsize_hires     ; log-10 spacing
    restwave = restwave_hires
    npix_rest = n_elements(restwave)

; build an observed-frame wavelength vector with constant 0.5-A
; spacing
    minwave_obs = 1250D
    maxwave_obs = 4D4
    dwave_obs = 0.5D
    npix_obs = (maxwave_obs-minwave_obs)/dwave_obs+1
    obswave = dindgen(npix_obs)*dwave_obs+minwave_obs

; build a sample of objects from the full AGES dataset that are
; suitable for template construction
    phot = read_ages(/photo)
    allzcat = read_ages(/ppxf)
    zcat = im_empty_structure(allzcat[0],ncopies=n_elements(phot),empty_value=-999.0)
    zcat[allzcat.ages_id] = allzcat
    isedfit = read_isedfit(isedfit_paramfile,thissfhgrid=sfhgrid,isedfit_dir=isedfit_dir)

    rejplates = [$
      104,$                     ; average fluxing vector
      106,110,209,310,311,$     ; unfluxed
      604]                      ; very crummy
    
;   splog, 'Hack!!'
;   index = where(zcat.pass eq 101 and zcat.continuum_snr gt 5 and $
;     zcat.z ge 0.05 and zcat.z le 0.85 and phot.imain,ngal)
    index = where(zcat.z ge 0.05 and zcat.z le 0.85 and phot.imain and $
      isedfit.chi2 lt 40 and phot.i_tot lt 19.95 and $
      (phot.pass ne rejplates[0]) and (phot.pass ne rejplates[1]) and $
      (phot.pass ne rejplates[2]) and (phot.pass ne rejplates[3]) and $
      (phot.pass ne rejplates[4]) and (phot.pass ne rejplates[5]) and $
      (phot.pass ne rejplates[6]),ngal)
;   index = index[0:20] & ngal = n_elements(index)
    splog, 'Total number of templates', ngal

; initialize the output data structure and header; the metadata header
; is written out below
    outinfo_rest = {$
      templateid:       0L,$
      d4000:           0.0,$
      hbeta_continuum: 0.0}
    outinfo_rest = replicate(outinfo_rest,ngal)
    outinfo_rest.templateid = lindgen(ngal)

    outinfo_obs = {$
      templateid:           0L,$
      pass:                  0,$
      aper:                  0,$
      ra:                   0D,$
      dec:                  0D,$
      z:                   0.0,$
      weight:              0.0,$
      ages_infiber_r:      0.0,$
      infiber_r:           0.0,$
      infiber_i:           0.0,$
      vdisp:               0.0,$ ; intrinsic velocity dispersion [km/s]
      sigma_kms:           0.0,$ ; intrinsic velocity linewidth [km/s]
      hbeta:               0.0,$
      hbeta_ew:            0.0,$
      oiii_hbeta:       -999.0,$
      decam_g:             0.0,$
      decam_r:             0.0,$
      decam_z:             0.0,$
      sb04_decam_r:        0.0,$ ; surface brightness [4" diameter aperture]
      logmstar:            0.0,$
      logsfr:              0.0,$
      av_ism:              0.0,$
      d4000:               0.0}
    outinfo_obs = replicate(outinfo_obs,ngal)
    outinfo_obs.templateid = lindgen(ngal)
    outinfo_obs.pass = zcat[index].pass
    outinfo_obs.aper = zcat[index].aper
    outinfo_obs.ra = phot[index].ra
    outinfo_obs.dec = phot[index].dec
    outinfo_obs.z = zcat[index].z ; heliocentric-corrected

    outinfo_obs.infiber_r = phot[index].infiber_r
    outinfo_obs.infiber_i = phot[index].infiber_i

    outinfo_obs.logmstar = isedfit[index].mstar_avg
    outinfo_obs.logsfr = isedfit[index].sfr_avg
    outinfo_obs.av_ism = isedfit[index].av*isedfit[index].mu ; mean attenuation of the ISM

; -------------------------
; get the galaxy weight (see build_mz_parent); in addition to the
; corrections for spectroscopic incompleteness (~2.1%), fiber
; incompleteness (~4.3%), and sparse-sampling, we also need to account
; for the catalog incompleteness (catalog_weight) and for the plates
; that were rejected above (field_weight); from Eisenstein, the
; catalog incompleteness is ~3%-6%, so assume an average value of 4%;
; finally, compute the final galaxy weight as the product of all the
; various selection terms                                                                           
    catalog_weight = 1.04
    allfield_weight = fltarr(ngal)+1

    allfield = fix(strmid(strtrim(outinfo_obs.pass,2),1,2)) ; field number
    field_weight = ages_upweight_field(rejplates,field=field)

    for ii = 0, n_elements(field)-1 do begin
       these = where(field[ii] eq allfield,nthese)
       if (nthese ne 0) then allfield_weight[these] = 1.0/field_weight[ii] ; note 1/FIELD_WEIGHT!
    endfor
    outinfo_obs.weight = phot[index].spec_weight*phot[index].target_weight*$
      phot[index].fiber_weight*catalog_weight*allfield_weight

; -------------------------
;; store the GALFIT parameters    
;    w1 = where(outinfo_obs.radius_halflight eq -999.0 and photo[index].flag_galfit_hi eq 0)
;    outinfo_obs[w1].radius_halflight = photo[index[w1]].re_galfit_hi*0.03 ; [arcsec]
;    outinfo_obs[w1].sersicn = photo[index[w1]].n_galfit_hi
;    outinfo_obs[w1].axis_ratio = photo[index[w1]].ba_galfit_hi
;
;    w2 = where(outinfo_obs.radius_halflight eq -999.0 and photo[index].flag_galfit_low eq 0)
;    outinfo_obs[w2].radius_halflight = photo[index[w2]].re_galfit_low*0.03 ; [arcsec]
;    outinfo_obs[w2].sersicn = photo[index[w2]].n_galfit_low
;    outinfo_obs[w2].axis_ratio = photo[index[w2]].ba_galfit_low

; -------------------------
; store the intrinsic velocity dispersion; if the velocity dispersion
; could not be modeled in AGES/PPXF then it was set to 165 km/s; fix
; that here 
;   ww = where(zcat.vdisp gt -900 and zcat.vdisp ne 165)
;   im_plothist, zcat[ww].vdisp, xx, yy, /peak, ysty=3
;   yfit = mpfitpeak(alog10(xx),yy,aa,/gauss,/pos,nterms=3)
;   djs_oplot, xx, yfit, color='green'

    aa = [0.94202718,2.1455627,0.18173899] ; =[peak,center,sigma] (log10)
    outinfo_obs.vdisp = zcat[index].vdisp
    fix = where(outinfo_obs.vdisp eq 165.0,nfix)
    if nfix ne 0L then outinfo_obs[fix].vdisp = (10.0^(randomn(seed,nfix)*(aa[2]-0.02)+aa[1]+0.05))<350.0
;   im_plothist, outinfo_obs.vdisp, bin=5, /peak
;   im_plothist, outinfo_obs[fix].vdisp, bin=5, color='red', /overplot, /peak

; -------------------------
; store the intrinsic emission-line width; if the velocity dispersion
; could not be modeled in AGES/PPXF then it was set to 0 km/s; fix
; that here 
;   ww = where(zcat.sigma_balmer gt -1.0)
;   im_plothist, zcat[ww].sigma_balmer, xx, yy, /peak, ysty=3, bin=5
;   yfit = mpfitpeak(alog10(xx),yy,aa,/gauss,/pos,nterms=3)
;   djs_oplot, xx, yfit, color='green', psym=10
    
    aa = [0.92850036,2.0729011,0.092070225] ; =[peak,center,sigma] (log10)
    outinfo_obs.sigma_kms = zcat[index].sigma_balmer
    fix = where(outinfo_obs.sigma_kms eq -1.0,nfix)
    if nfix ne 0L then outinfo_obs[fix].sigma_kms = (10.0^(randomn(seed,nfix)*(aa[2])+aa[1]+0.04))<400.0

    broad = where(zcat[index].isbroad)
    outinfo_obs[broad].sigma_kms = zcat[index[broad]].sigma_broad
    
;   im_plothist, outinfo_obs.sigma_kms, bin=5, /peak
;   im_plothist, outinfo_obs[fix].sigma_kms, bin=5, color='red', /overplot, /peak

; -------------------------
; store the H-beta emission-line EW
    good = where(zcat[index].h_beta_ew[1] gt 0,ngood,comp=limit,ncomp=nlimit)
    if ngood ne 0L then outinfo_obs[good].hbeta_ew = zcat[index[good]].h_beta_ew[0]
    if nlimit ne 0L then outinfo_obs[limit].hbeta_ew = zcat[index[limit]].h_beta_ew_limit

; -------------------------
; store the [OIII]/H-beta ratio
;   ww = where(zcat[index].oiii_5007[1] gt 0 and zcat[index].h_beta[1] gt 0)
;   oiiihb = alog10(zcat[index[ww]].oiii_5007[0]/zcat[index[ww]].h_beta[0])
;   im_plothist, oiiihb, xx, yy, /peak, ysty=3, bin=0.05
;   yfit = mpfitpeak(xx,yy,aa,/gauss,/pos,nterms=3)
;   djs_oplot, xx, yfit, color='green', psym=10

    aa = [0.96964753,-0.13746193,0.33575511]
    good = where(zcat[index].oiii_5007[1] gt 0 and zcat[index].h_beta[1] gt 0,ngood)
    outinfo_obs[good].oiii_hbeta = alog10(zcat[index[good]].oiii_5007[0]/zcat[index[good]].h_beta[0])

    fix = where(outinfo_obs.oiii_hbeta eq -999.0,nfix)
    if nfix ne 0L then outinfo_obs[fix].oiii_hbeta = randomn(seed,nfix)*aa[2]+aa[1]

;   im_plothist, outinfo_obs.oiii_hbeta, bin=0.05, /peak
;   im_plothist, outinfo_obs[fix].oiii_hbeta, bin=0.05, color='red', /overplot, /peak
    
; initialize the output flux vectors
    outflux_obs = fltarr(ngal,npix_obs)
    outflux_rest = fltarr(ngal,npix_rest)
    outflux_continuum_rest = fltarr(ngal,npix_rest)

; sort by CHUNKINDX for speed since all the overhead is in
; READ_ISEDFIT() 
    allchunkindx = isedfit[index].chunkindx
    chunkindx = allchunkindx[uniq(allchunkindx,sort(allchunkindx))]
    nchunk = n_elements(chunkindx)
    
;   for ichunk = 1, nchunk-1 do begin
    for ichunk = 0, nchunk-1 do begin
       splog, format='("Working on chunk ",I0,"/",I0)', ichunk+1, nchunk
       these = where(chunkindx[ichunk] eq allchunkindx)
       if keyword_set(debug) then these = these[0:10]
       nthese = n_elements(these)

; (re)build the iSEDfit models
       ised = read_isedfit(isedfit_paramfile,/getmodels,isedfit_dir=isedfit_dir,$
         montegrids_dir=isedfit_dir+'montegrids/',index=index[these],/flam,$
         noigm=0,thissfhgrid=sfhgrid) 

; read the AGES 1D spectra and best-fitting models       
       spec1d = read_ages_gandalf_specfit(zcat[index[these]],$
         /solar,/observed,/silent)

;; rebuild the model spectrum over a wider wavelength range       
;       modelflux = restore_ages_ppxf_bestfit(zcat[index[these]].continuum_coeff,$
;         ebv=zcat[index[these]].continuum_ebv,bestwave=modelwave,bestage=modelage,$
;         /solar,inst_vdisp=ages_ppxf_instvdisp())
       
;      for igal = 20, 20 do begin
       for igal = 0, nthese-1 do begin
          print, format='("Building template ",I0,"/",I0,A10,$)', igal, nthese, string(13b)

          zobj = outinfo_obs[these[igal]].z
          dlum = pc10*10D^(lf_distmod(zobj,omega0=0.3D,omegal0=0.7D)/5D)/0.7D
          distratio = (dlum/pc10)^2 ; distance ratio factor (obs-->rest)

; calculate the fraction of the flux within the fiber          
          outinfo_obs[these[igal]].ages_infiber_r = ((k_project_filters(k_lambda_to_edges($
            exp(spec1d[igal].wave)),spec1d[igal].flux,filterlist='decam_r.par')/$
            k_project_filters(k_lambda_to_edges(ised[igal].wave),ised[igal].flux,$
            filterlist='decam_r.par'))<1.0)>0.05 ; note!

; get the rest-frame *continuum* flux around H-beta and thereby the
; integrated emission-line flux
          cflux_hbeta = isedfit_linecontinuum(ised[igal].wave/(1+zobj),$
            ised[igal].flux*(1+zobj),linewave=4861.32,debug=0) ; wavelength in air
          outinfo_rest[these[igal]].hbeta_continuum = cflux_hbeta
          outinfo_obs[these[igal]].hbeta = cflux_hbeta*outinfo_obs[these[igal]].hbeta_ew ; [erg/s/cm2]

; measure D(4000) from the emission-line free (rest-frame) continuum
; spectrum 
          outinfo_obs[these[igal]].d4000 = im_d4000(ised[igal].wave/(1+zobj),$
            (ised[igal].flux-ised[igal].nebflux)*(1+zobj),debug=0)
          outinfo_rest[these[igal]].d4000 = outinfo_obs[these[igal]].d4000

; construct a rest-frame emission-line spectrum which has the correct
; H-beta flux; call ISEDFIT_NEBULAR just so we can get the LINE
; structure 
          linesigma = outinfo_obs[these[igal]].sigma_kms
          reflineflux = outinfo_obs[these[igal]].hbeta ; [erg/s/cm2]

          junk = isedfit_nebular(10D^ised[igal].nlyc,inst_vsigma=0D,vsigma=linesigma,$
            oiiihb=10D^outinfo_obs[these[igal]].oiii_hbeta,wave=10D^restwave,line=line,$
            oiidoubletratio=oiidoubletratio,/vacuum,/nospectrum)

; rescale all the lines to the H-beta flux @10pc, normalized to 1Msun;
; note that there is no (1+z) factor because this is an integrated
; flux
          ishbeta = where(strtrim(line.name,2) eq 'Hbeta')
;         hbetafactor = reflineflux/line[ishbeta].flux
          hbetafactor = reflineflux/line[ishbeta].flux*distratio/10D^ised[igal].mstar ; no (1+z)!
          for ll = 0, n_elements(line)-1 do begin
             line[ll].flux = line[ll].flux*hbetafactor
             line[ll].amp = line[ll].amp*hbetafactor
          endfor

;; test to be sure we get the right total H-beta flux - yes!
;          ishbeta = where(strtrim(line.name,2) eq 'Hbeta')
;          testemspectrum = build_emline(restwave*0,logwave=restwave,$
;            lineflux=line[ishbeta].flux,linesigma=linesigma,linewave=line[ishbeta].wave)
;          print, im_integral(10^restwave,testemspectrum)*10D^ised[igal].mstar/distratio, $
;            outinfo_obs[these[igal]].hbeta, linesigma

; now build the full emission-line spectrum at the new velocity
; resolution
          emspectrum_rest = restwave*0
          for ll = 0, n_elements(line)-1 do emspectrum_rest = build_emline(emspectrum_rest,$
            logwave=restwave,lineflux=line[ll].flux,linesigma=linesigma,linewave=line[ll].wave)

; construct the continuum rest-frame spectrum (@10pc, normalized to 1
; Msun), being sure to subtract out the emission-line spectrum used by
; iSEDfit; also convolve the model to the inferred stellar velocity
; dispersion
          airtovac, ised[igal].wave/(1+zobj), isedwave_rest ; [Angstrom]
          isedflux_rest = (ised[igal].flux-ised[igal].nebflux)*(1+zobj)*$ ; [erg/s/cm2/A]
            distratio/10D^ised[igal].mstar

          continuum_rest1 = interpol(isedflux_rest,alog10(isedwave_rest),restwave)
          continuum_rest = continuum_convolve(continuum_rest1,velscale=velpixsize_hires,$
            vdisp=outinfo_obs[these[igal]].vdisp)

          restflux = continuum_rest+emspectrum_rest

; now build the observed-frame spectrum
          airtovac, ised[igal].wave, isedwave_obs           ; [Angstrom]
          isedflux_obs = ised[igal].flux-ised[igal].nebflux ; [erg/s/cm2/A]

          continuum_obs1 = interpol(isedflux_obs,isedwave_obs,obswave)
          continuum_obs = continuum_convolve(continuum_obs1,velscale=velpixsize_hires,$
            vdisp=outinfo_obs[these[igal]].vdisp)

;         emspectrum_obs = interpol(emspectrum_rest/(1+zobj),10D^restwave*(1+zobj),obswave)
          emspectrum_obs = interpol(emspectrum_rest*10D^ised[igal].mstar/(1+zobj)/distratio,$
            10D^restwave*(1+zobj),obswave)
          obsflux = continuum_obs+emspectrum_obs
          
; pack it in!
          outflux_rest[these[igal],*] = restflux
          outflux_continuum_rest[these[igal],*] = continuum_rest
          outflux_obs[these[igal],*] = obsflux
;         outflux.flux[igal,*] = flux

; synthesize photometry and pack in the output structure
          obswave_edges = k_lambda_to_edges(obswave)
          grz = k_project_filters(obswave_edges,obsflux,filterlist=synth_filters)
          outinfo_obs[these[igal]].decam_g = -2.5*alog10(reform(grz[0,0,0]))
          outinfo_obs[these[igal]].decam_r = -2.5*alog10(reform(grz[0,0,1]))
          outinfo_obs[these[igal]].decam_z = -2.5*alog10(reform(grz[0,0,2]))
;         if igal eq 50 then stop

; get the r-band surface brightness
          outinfo_obs[these[igal]].sb04_decam_r = phot[index[these[igal]]].r_mag_aper_04 + $
            (outinfo_obs[these[igal]].decam_r+2.5*alog10(reform(k_project_filters($
            obswave_edges,obsflux,filterlist='ndwfs_R.par')))) + 2.5*alog10(2.0*!pi*4.0^2)

; debugging plot          
          if keyword_set(debug) then begin
             xrange = [1200,6E4]
             inrange = where(ised[igal].wave gt xrange[0] and ised[igal].wave lt xrange[1])
             good = where(spec1d[igal].wave ne 0)
             yrange = minmax(ised[igal].flux[inrange])*[0.6,2]
;            yrange = [min(ised[igal].flux[inrange])*0.9,max(spec1d[igal].flux[good])*1.1]
             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
               xrange=xrange/1D4, yrange=yrange, /xlog, /ylog
             djs_oplot, exp(spec1d[igal].wave[good])/1D4, spec1d[igal].flux[good]
;            djs_oplot, modelwave*(1+zobj), modelflux[*,igal]/(1+zobj), color=cgcolor('dodger blue')
             djs_oplot, exp(spec1d[igal].wave[good])/1D4, spec1d[igal].continuum[good]+$
               spec1d[igal].linefit[good], color='red'
             djs_oplot, ised[igal].wave/1D4, ised[igal].flux, color='blue'
;            djs_oplot, ised[igal].wave/1D4, ised[igal].flux*outinfo_obs[these[igal]].infiber_r, $
;              color=cgcolor('dodger blue')
             djs_oplot, obswave/1D4, obsflux*outinfo_obs[these[igal]].infiber_r, $
               color=cgcolor('forest green')

             djs_oplot, weff/1D4, isedfit[index[these[igal]]].maggies*maggies2flam, $
               psym=6, symsize=4, color='yellow'
             im_legend, [string(zcat[index[these[igal]]].pass,format='(I3.3)')+'/'+$
               string(zcat[index[these[igal]]].aper,format='(I3.3)'),'z='+$
               strtrim(string(zcat[index[these[igal]]].z,format='(F12.4)'),2)], $
               /left, /top, box=0
; inset
             xrange = [4000,6000]
             inrange = where(ised[igal].wave gt xrange[0] and ised[igal].wave lt xrange[1])
             yrange = [min(ised[igal].flux[inrange])<im_min(spec1d[igal].flux[good],sigrej=2),$
               max(ised[igal].flux[inrange])>im_max(spec1d[igal].flux[good],sigrej=2)]*[0.9,1.1]
             djs_plot, [0], [0], /nodata, /noerase, position=[0.4,0.15,0.8,0.5], $
               /norm, xsty=1, ysty=1, xrange=xrange/1D4, yrange=yrange, $
               ytickname=replicate(' ',10), charsize=1.3
             djs_oplot, exp(spec1d[igal].wave[good])/1D4, spec1d[igal].flux[good]
;            djs_oplot, ised[igal].wave/1D4, ised[igal].flux*outinfo_obs[these[igal]].infiber_r, $
;              color=cgcolor('dodger blue')
             djs_oplot, exp(spec1d[igal].wave[good])/1D4, spec1d[igal].continuum[good]+$
               spec1d[igal].linefit[good], color='red'
             djs_oplot, obswave/1D4, obsflux*outinfo_obs[these[igal]].infiber_r, color=cgcolor('forest green')
             cc = get_kbrd(1)
          endif
       endfor                   ; close galaxy loop
    endfor                      ; close chunk loop

; initialize the rest-frame and observed-frame spectral headers
    mkhdr, outhdr_rest, transpose(outflux_rest), /extend
    sxdelpar, outhdr_rest, 'DATE'
    sxdelpar, outhdr_rest, 'COMMENT'
    sxaddpar, outhdr_rest, 'VERSION', version, ' template version number'
    sxaddpar, outhdr_rest, 'OBJTYPE', 'BGS-AGES', ' object type'
    sxaddpar, outhdr_rest, 'DISPAXIS', 1, ' dispersion axis'
    sxaddpar, outhdr_rest, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
    sxaddpar, outhdr_rest, 'CUNIT1', 'Angstrom', ' units of wavelength array'
    sxaddpar, outhdr_rest, 'CRPIX1', 1, ' reference pixel number'
    sxaddpar, outhdr_rest, 'CRVAL1', min(restwave), ' reference log10(Angstrom)'
    sxaddpar, outhdr_rest, 'CDELT1', pixsize_hires, ' delta log10(Angstrom)'
    sxaddpar, outhdr_rest, 'LOGLAM', 1, ' log10 spaced wavelengths?'
    sxaddpar, outhdr_rest, 'WAVEMIN', min(10D^restwave), ' minimum wavelength (Angstrom)'
    sxaddpar, outhdr_rest, 'WAVEMAX', max(10D^restwave), ' maximum wavelength (Angstrom)'
    sxaddpar, outhdr_rest, 'WAVEUNIT', 'Angstrom', ' wavelength units'
    sxaddpar, outhdr_rest, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
    sxaddpar, outhdr_rest, 'VELSCALE', velpixsize_hires, ' pixel size in km/s'
    sxaddpar, outhdr_rest, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    sxaddpar, outhdr_rest, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'

    mkhdr, outhdr_continuum_rest, transpose(outflux_continuum_rest), /extend
    sxdelpar, outhdr_continuum_rest, 'DATE'
    sxdelpar, outhdr_continuum_rest, 'COMMENT'
    sxaddpar, outhdr_continuum_rest, 'VERSION', version, ' template version number'
    sxaddpar, outhdr_continuum_rest, 'OBJTYPE', 'BGS-AGES', ' object type'
    sxaddpar, outhdr_continuum_rest, 'DISPAXIS', 1, ' dispersion axis'
    sxaddpar, outhdr_continuum_rest, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
    sxaddpar, outhdr_continuum_rest, 'CUNIT1', 'Angstrom', ' units of wavelength array'
    sxaddpar, outhdr_continuum_rest, 'CRPIX1', 1, ' reference pixel number'
    sxaddpar, outhdr_continuum_rest, 'CRVAL1', min(restwave), ' reference log10(Angstrom)'
    sxaddpar, outhdr_continuum_rest, 'CDELT1', pixsize_hires, ' delta log10(Angstrom)'
    sxaddpar, outhdr_continuum_rest, 'LOGLAM', 1, ' log10 spaced wavelengths?'
    sxaddpar, outhdr_continuum_rest, 'WAVEMIN', min(10D^restwave), ' minimum wavelength (Angstrom)'
    sxaddpar, outhdr_continuum_rest, 'WAVEMAX', max(10D^restwave), ' maximum wavelength (Angstrom)'
    sxaddpar, outhdr_continuum_rest, 'WAVEUNIT', 'Angstrom', ' wavelength units'
    sxaddpar, outhdr_continuum_rest, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
    sxaddpar, outhdr_continuum_rest, 'VELSCALE', velpixsize_hires, ' pixel size in km/s'
    sxaddpar, outhdr_continuum_rest, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    sxaddpar, outhdr_continuum_rest, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'

    mkhdr, outhdr_obs, transpose(outflux_obs), /extend
    sxdelpar, outhdr_obs, 'DATE'
    sxdelpar, outhdr_obs, 'COMMENT'
    sxaddpar, outhdr_obs, 'VERSION', version, ' template version number'
    sxaddpar, outhdr_obs, 'OBJTYPE', 'BGS-AGES', ' object type'
    sxaddpar, outhdr_obs, 'DISPAXIS', 1, ' dispersion axis'
    sxaddpar, outhdr_obs, 'CTYPE1', 'WAVE-WAV', ' meaning of wavelength array'
    sxaddpar, outhdr_obs, 'CUNIT1', 'Angstrom', ' units of wavelength array'
    sxaddpar, outhdr_obs, 'CRPIX1', 1, ' reference pixel number'
    sxaddpar, outhdr_obs, 'CRVAL1', min(obswave), ' reference wavelength (Angstrom)'
    sxaddpar, outhdr_obs, 'CDELT1', dwave_obs, ' dispersion (Angstrom)'
    sxaddpar, outhdr_obs, 'LOGLAM', 0, ' log10 spaced wavelengths?'
    sxaddpar, outhdr_obs, 'WAVEMIN', min(obswave), ' minimum wavelength (Angstrom)'
    sxaddpar, outhdr_obs, 'WAVEMAX', max(obswave), ' maximum wavelength (Angstrom)'
    sxaddpar, outhdr_obs, 'WAVEUNIT', 'Angstrom', ' wavelength units'
    sxaddpar, outhdr_obs, 'AIRORVAC', 'vac', ' wavelengths in vacuum (vac) or air'
    sxaddpar, outhdr_obs, 'BUNIT', 'erg/s/cm2/A', ' spectrum flux units'
    sxaddpar, outhdr_obs, 'FLUXUNIT', 'erg/s/cm2/A', ' spectrum flux units'

; we need to do some MWRFITS jujitsu to get the metadata header right 
    units_obs = [$
      'TEMPLATEID,,unique template ID number (0-indexed)',$
      'PASS,,AGES pass number',$
      'APER,,AGES aperture number (0-indexed)',$
      'RA,degrees,right ascension (J2000)',$
      'DEC,degrees,declination (J2000)',$
      'Z,,AGES heliocentric redshift',$
      'WEIGHT,,statistical weight',$
      'INFIBER_R,,fraction of r-band light in the AGES fiber',$
      'VDISP,km/s,stellar velocity dispersion',$
      'SIGMA_KMS,km/s,emission line velocity width',$
      'HBETA,erg/s/cm2,H-beta emission-line flux',$
      'HBETA_EW,Angstrom,rest-frame H-beta emission-line equivalent width',$
      'OIII_HBETA,,logarithmic [OIII] 5007/H-beta flux ratio',$
      'DECAM_G,mag,synthesized DECam g-band AB mag',$
      'DECAM_R,mag,synthesized DECam r-band AB mag',$
      'DECAM_Z,mag,synthesized DECam z-band AB mag',$
      'SB04_DECAM_R,mag/arcsec2,r-band surface brightness (4 arcsec diameter)',$
      'LOGMSTAR,Msun,log10(stellar mass) (Chabrier, h=0.7)',$
      'LOGSFR,Msun/yr,log10(star formation rate) (Chabrier, h=0.7)',$
      'AV_ISM,mag,V-band attenuation (Charlot & Fall 2000)',$
      'D4000,,4000-Angstrom break']

    units_rest = [$
      'TEMPLATEID,,unique template ID number (0-indexed)',$
      'D4000,,4000-Angstrom break',$
      'HBETA_CONTINUUM,erg/s/cm2/A,rest-frame continuum flux around H-beta']
    units_continuum_rest = units_rest

    newheader = [$
;     'BGS/AGES galaxy templates',$
;     'Generated by J. Moustakas (jmoustakas@siena.edu) on '+im_today(),$
      'EXTNAME'+" = '"+string('METADATA','(A-8)')+"'           /"]

; observed frame
    newheader_obs = newheader
    mwrfits, 0, outfile_obs, /create
    mwrfits, outinfo_obs, outfile_obs, /silent
    metahdr_obs =  im_update_header(outfile_obs,units_obs,newheader_obs)
    
    mwrfits, transpose(outflux_obs), outfile_obs, outhdr_obs, /create
    mwrfits, outinfo_obs, outfile_obs, metahdr_obs

; rest-frame - full spectrum
    newheader_rest = newheader
    mwrfits, 0, outfile_rest, /create
    mwrfits, outinfo_rest, outfile_rest, /silent
    metahdr_rest =  im_update_header(outfile_rest,units_rest,newheader_rest)

    im_mwrfits, transpose(outflux_rest), outfile_rest, outhdr_rest, /clobber, /nogzip
    im_mwrfits, outinfo_rest, outfile_rest, metahdr_rest, /append, /nogzip

; rest-frame - continuum spectrum
    newheader_continuum_rest = newheader
    mwrfits, 0, outfile_continuum_rest, /create
    mwrfits, outinfo_continuum_rest, outfile_continuum_rest, /silent
    metahdr_continuum_rest =  im_update_header(outfile_continuum_rest,$
      units_continuum_rest,newheader_continuum_rest)

    im_mwrfits, transpose(outflux_continuum_rest), outfile_continuum_rest, $
      outhdr_continuum_rest, /clobber, /nogzip
    im_mwrfits, outinfo_rest, outfile_continuum_rest, metahdr_continuum_rest, /append, /nogzip

stop    
    
return
end
