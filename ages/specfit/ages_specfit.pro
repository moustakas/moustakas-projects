;+
; NAME:
;   AGES_SPECFIT
;
; PURPOSE:
;   Fit the fluxed and unfluxed AGES spectra.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   thispass  - explicit pass name to fit (e.g., 'ages_101.fits.gz')
;   firstpast - first pass number to fit 
;   lastpass  - last pass number to fit
;   index     - index number of object to fit (for testing)
;
; KEYWORD PARAMETERS: 
;   unfluxed - fit the unfluxed plates
;   notweak  - do not use the 'tweaked' AGES spectra
;   solar    - only use the solar-metallicity templates (default is to
;     use three metallicities unless /NOTWEAK or /UNFLUXED)
;   test     - generate 'test' files
;   doplot   - do not generate output; render the plots on the screen 
;   refitbad - refit spectra where the line-fits were 'bad'; these are
;     objects where chi2>5 on any of [OII], H-beta, [OIII], H-alpha,
;     or [NII]; at the same time, relax the ZINDEX constraints to
;     allow for more freedom, and force /SOLAR
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2004 Jun 04, U of A - originally written 
;   jm04sep12uofa - updated to fit v0300; also we now fit each plate 
;     separately since the fluxing issues vary across plates 
;   jm05jul30uofa - updated for release 2.0
;   jm06feb05uofa - fit each plate separately; fit the post-skysubpca 
;                   galaxy spectra, including the QSO's
;   jm08apr03nyu - optionally do not repair the spectra for poor
;    fluxing  
;   jm08aug11nyu - a bit of a major rewrite 
;   jm08oct27nyu - added REFITBAD keyword
;
; Copyright (C) 2004-2006, 2008, John Moustakas
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

pro ages_specfit, thispass=thispass, firstpass=firstpass, lastpass=lastpass, $
  unfluxed=unfluxed, notweak=notweak, solar=solar, test=test, doplot=doplot, $
  index=index1, refitbad=refitbad

    if keyword_set(doplot) or keyword_set(test) then $
      silent = 0 else silent = 1

; ssh -X prism
; echo "ages_specfit, /notweak" | idl > & ages_specfit.log.1 &

    base_specfitpath = ages_path(/specfit)
    base_spec1dpath = ages_path(/spec1d)

    if keyword_set(unfluxed) then begin
       version = ages_version(/unfluxed_specfit)
       specfitpath = base_specfitpath+'unfluxed/'+version+'/'
       spec1dpath = base_spec1dpath+'unfluxed/after_skysubpca/'
       linefile = base_specfitpath+'elinelist_unfluxed_'+version+'.dat'
    endif else begin
       version = ages_version(/ispec_specfit)
       if keyword_set(notweak) then begin
          specfitpath = base_specfitpath+'fluxed/notweak/'+version+'/'
          spec1dpath = base_spec1dpath+'fluxed/after_skysubpca/' 
       endif else begin
          specfitpath = base_specfitpath+'fluxed/tweak/'+version+'/'
          spec1dpath = base_spec1dpath+'fluxed/tweak/'
       endelse
       if keyword_set(refitbad) then $
         linefile = base_specfitpath+'elinelist_untied_'+version+'.dat' else $
           linefile = base_specfitpath+'elinelist_'+version+'.dat'
    endelse

    indexfile = base_specfitpath+'indexlist_'+version+'.dat'
    if keyword_set(refitbad) then solar = 1 ; NOTE!
    if keyword_set(notweak) or keyword_set(solar) then $
      Zbc03 = ['Z02'] else $
      Zbc03 = ['Z004','Z02','Z05']
    
    templatefile = base_specfitpath+'BC03_'+Zbc03+'_salpeter_ages_'+$
      ages_version(/templates)+'_templates.fits'
       
    if (n_elements(thispass) eq 0L) then $
      thispass = file_search(spec1dpath+'ages_???.fits.gz') else $
      thispass = spec1dpath+thispass
    npass = n_elements(thispass)

    if (n_elements(firstpass) eq 0L) then firstpass = 0L
    if (n_elements(lastpass) eq 0L) then lastpass = npass-1L
    
; fit each plate separately
    t0 = systime(1)
    for ipass = firstpass, lastpass do begin

; read the spectra
       splog, 'Reading '+thispass[ipass]
       specpass = mrdfits(thispass[ipass],1,/silent)
       pass1 = specpass.pass[0]
       suffix = string(pass1,format='(I0)')       
       if keyword_set(refitbad) then suffix = suffix+'_refitbad'
       if keyword_set(test) then suffix = suffix+'_test'

       npix = n_elements(specpass.wave)
       nobj = n_elements(specpass.z)
       
       prepend = {ages_id: -1L, galaxy: '', ra: 0.0D, $
         dec: 0.0D, pass: -1, aper: -1, class: '', $
         z: 0.0, spzbest: 0.0, fluxing_problem: -1, $
         zmerge_multiple: -1}
       prepend = replicate(prepend,nobj)
       prepend.ages_id         = specpass.ages_id
       prepend.galaxy          = specpass.galaxy
       prepend.ra              = specpass.ra
       prepend.dec             = specpass.dec
       prepend.pass            = specpass.pass
       prepend.aper            = specpass.aper
       prepend.class           = specpass.class
       prepend.z               = specpass.z
       prepend.spzbest         = specpass.spzbest
       prepend.fluxing_problem = specpass.fluxing_problem
       prepend.zmerge_multiple = specpass.zmerge_multiple

; build an inverse variance cube for ispeclinefit(); note that double
; precision is crucial, especially for invvar!

       galaxy = strtrim(specpass.galaxy,2)
       class = strtrim(specpass.class,2)
       zobj = specpass.z
       flux = specpass.flux
       wave = rebin(reform(specpass.wave,npix,1),npix,nobj)
       invvar = (1.0D/(specpass.ferr+(specpass.ferr eq 0.0))^2.0)*$
         (specpass.ferr ne 0.0)

;      check = where((finite(invvar) eq 0) or (finite(specpass.ferr) eq 0),ncheck)
;      if (ncheck ne 0L) then message, 'Uh-oh!'
;      continue

       if keyword_set(refitbad) then begin
          prevfitfile = file_search(specfitpath+'?????_'+string(pass1,format='(I0)')+$
            '_specdata.fits.gz',count=nprevcount)
          if (nprevcount ne 1L) then message, 'Handle this case'
          splog, 'Reading '+prevfitfile
          prev = mrdfits(prevfitfile,1)
          case pass1 of
             214: begin
                prevbad = speclinefit_locate(prev,'ages_214/108') ; requires /SOLAR ; pg 108
                vlinemaxtweak = 300.0 ; required by 214/108
                if (Zbc03 ne 'Z02') then message, 'Fitting 214/108 requires solar metallicity!'
             end
             else: begin
                delvarx, vlinemaxtweak
                splog, 'No refitting necessary for pass '+string(pass1,format='(I0)')+'!'
                return
             end
          endcase
          nprevbad = n_elements(prevbad)
;         prevbad = where((prev.z gt 0.01) and (prev.z lt 1.0) and $
;           (((prev.oii_3727[1] gt 0.0) and (prev.oii_3727_chi2 gt 5.0)) or $
;           ((prev.oiii_5007[1] gt 0.0) and (prev.oiii_5007_chi2 gt 5.0)) or $
;           ((prev.nii_6584[1] gt 0.0) and (prev.nii_6584_chi2 gt 5.0)) or $
;           ((prev.h_alpha[1] gt 0.0) and (prev.h_alpha_chi2 gt 5.0)) or $
;           ((prev.h_beta[1] gt 0.0) and (prev.h_beta_chi2 gt 5.0))),nprevbad)
          if (nprevbad ne 0L) then begin
             prepend = prepend[prevbad]
             wave = wave[*,prevbad]
             flux = flux[*,prevbad]
             invvar = invvar[*,prevbad]
; use the absorption-line redshift as the initial guess! so that
; VLINEMAXTWEAK can be a smaller quantity
;            zobj = prev[prevbad].z_abs 
             zobj = zobj[prevbad]
             galaxy = galaxy[prevbad]
             class = class[prevbad]
             nobj = nprevbad
          endif else begin
             splog, 'No bad objects to refit for pass, '+string(pass1,format='(I0)')+'!'
             continue
          endelse
       endif

       if (n_elements(index1) eq 0L) then index = lindgen(nobj) else $
         index = index1

       if keyword_set(unfluxed) then begin

;         index = where(strtrim(class,2) eq 'QSO')
          
          specres = 4.5 ; 5.0        ; FWHM instrumental resolution [A]
          vlinemaxtweak = 500.0 ; [km/s]
          sigmax = replicate(5000.0,nobj)
          
          specdata = ispeclinefit_unfluxed(wave[*,index],flux[*,index],$
            invvar[*,index],zobj=zobj[index],specres=specres,sigmax=sigmax[index],     $
            galaxy=galaxy[index],outpath=specfitpath,suffix=suffix,$
            linefile=linefile,specfit=specfit,/nologfile,vlinemaxtweak=vlinemaxtweak,$
            /clobber,silent=silent,doplot=doplot,specdatafile=specdatafile)
          
; prepend some other info and overwrite SPECDATAFILE

          mwrfits, struct_addtags(prepend[index],struct_trimtags($
            specdata,except=['galaxy'])), specdatafile, /create
          spawn, 'gzip -f '+specdatafile, /sh

       endif else begin
          
; specify some of the parameters of the fitting, which were determined
; by trial and error; the relatively low SPECRES is needed because
; some lines are narrower than the instrumental resolution; make up
; for this by upping the velocity dispersion
          
          continuum_maxwave = 8500.0 ; exclude wavelengths contaminated by the red leak
          specres = 4.5 ; 5.0        ; FWHM instrumental resolution [A]
          vdisp = 50.0               ; assume constant stellar velocity dispersion [km/s]

; think hard about making these different; in general, the emission
; lines should be able to "escape" back to the input guess redshift,
; since we define all the rest-frame quantities in ispeclinefit()
; according to z_abs          
          
          if (n_elements(vmaxtweak) eq 0L) then vmaxtweak = 300.0         ; [km/s]
          if (n_elements(vlinemaxtweak) eq 0L) then vlinemaxtweak = 500.0 ; 700 [km/s]

          sigmax = replicate(500.0,nobj)
          qso = where(strtrim(class,2) eq 'QSO',nqso)
          if (nqso ne 0L) then sigmax[qso] = 5000.0 ; [km/s]

;         zsnrcross = replicate(-1.0,nobj)
          zsnrcross = replicate(2.0,nobj)
          hiz = where(zobj gt 1.0,nhiz)
          if (nhiz ne 0L) then zsnrcross[hiz] = -1.0 ; suppress

; specify MINWAVE for plates that have problematic spectrophotometry
; in the blue so that the fits aren't affected; do this only for
; /NOTWEAK; note that 104 & 312 are the least well fluxed plates 

          if keyword_set(notweak) then begin
             case pass1 of
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
          endif
             
; now do the fit

          specdata = ispeclinefit(wave[*,index],flux[*,index],invvar[*,index],$
            zobj=zobj[index],specres=specres,vdisp=vdisp,sigmax=sigmax[index],$
            zsnrcross=zsnrcross[index],galaxy=galaxy[index],outpath=specfitpath,$
            suffix=suffix,templatefile=templatefile,linefile=linefile,indexfile=indexfile,$
            specfit=specfit,continuum_minwave=continuum_minwave,$
            continuum_maxwave=continuum_maxwave,/nologfile,vmaxtweak=vmaxtweak,$
            vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,doplot=doplot,$
            specdatafile=specdatafile)
             
; prepend some other info and overwrite SPECDATAFILE

          mwrfits, struct_addtags(prepend[index],struct_trimtags($
            specdata,except=['galaxy'])), specdatafile, /create
          spawn, 'gzip -f '+specdatafile, /sh

       endelse 

    endfor 
    splog, 'Total '+string((systime(1)-t0)/3600.0,format='(G0.0)')+' hours'

; merge all the plates and write out

stop    
    
    ages_merge_specfit, alldata, unfluxed=unfluxed, notweak=notweak, $
      test=test, /write
    
stop
    
return
end    
