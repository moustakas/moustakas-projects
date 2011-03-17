pro sings_specfit, sings1, drift56=drift56, drift20=drift20, $
  nuclear=nuclear, test=test, doplot=doplot, index=index1
; jm05jul26uofa
; jm08sep08nyu - major update begun
; jm08oct19nyu - major update completed

    if (n_elements(sings1) eq 0L) then sings1 = sings_read_info()

    version = sings_version(/specfit)
    spec1dpath = sings_path(/spec1d)
    base_specfitpath = sings_path(/specfit)
    specfitpath = base_specfitpath+version+'/'

    Zbc03 = ['Z02'] ; ['Z004','Z02','Z05']
    templatefile = base_specfitpath+'BC03_'+Zbc03+'_salpeter_sings_'+$
      sings_version(/templates)+'_templates.fits'

    linefile = base_specfitpath+'elinelist_'+version+'.dat'
    linefile_broad = base_specfitpath+'elinelist_broad_'+version+'.dat'
    indexfile = base_specfitpath+'indexlist_'+version+'.dat'

; initialize the broad-line structure

    linefit_broad = icreate_linefit(3)
    linefit_broad.linename = ['H_alpha_broad','H_beta_broad','H_gamma_broad']
    linefit_broad.linewave = [6562.80,4861.325,4340.464]
    linefit_broad = struct_trimtags(parse_ilinefit(linefit_broad),$
      except=['LINENAME','Z_LINE*'])

; ##################################################    
; fit the nuclear spectra

    if keyword_set(nuclear) then begin

       if keyword_set(test) then suffix = 'test_nuclear' else suffix = 'sings_nuclear'

       nuc = where(sings1.nuclear,nnuc)
       sings = sings1[nuc]
       if (n_elements(index1) ne 0L) then sings = sings1[nuc[index1]]
       nnuc = n_elements(sings)

       speclist = strtrim(sings.nuclear_file,2)
       galaxy = strtrim(sings.galaxy,2)

; read the spectra; linearly extrapolate the wavelength vector to
; avoid NaN in ibackfit 

       forage = iforage(spec1dpath+speclist)
       npix = max(forage.naxis1)
       wave = dblarr(npix,nnuc)
       flux = wave*0.0D
       invvar = wave*0.0D
       for inuc = 0L, nnuc-1L do begin
          spec1 = rd1dspec(speclist[inuc],datapath=spec1dpath)
          wave[*,inuc] = interpol(spec1.wave,lindgen(spec1.npix),lindgen(npix))
;         wave[0L:spec1.npix-1L,inuc] = spec1.wave
          flux[0L:spec1.npix-1L,inuc] = spec1.spec
          invvar[0L:spec1.npix-1L,inuc] = 1.0/spec1.sigspec^2.0
          icleanup, spec1
       endfor
       
; note that values of VMAXTWEAK smaller than about 1000 km/s are
; useless because the pixels are themselves about ~160 km/s big
       
; read in the parameters of the fitting       

       nucparams1 = rsex(specfitpath+'zguess_nuclear.sex')
       these = speclinefit_locate(sings,nucparams1.galaxy)
       good = where(these ne -1,ngood)
       if (ngood ne 0L) then nucparams = nucparams1[good]

       vdisp = nucparams.vdisp
       vmaxtweak = nucparams.vmaxtweak
       vlinemaxtweak = nucparams.vlinemaxtweak
       sigmax = nucparams.sigmax
       zsnrcross = nucparams.zsnrcross
       zobj = nucparams.zobj
       
       specres = 7.5 ; [FWHM, A]

; now fit them!       

       specdata_all = ispeclinefit(wave,flux,invvar,zobj=zobj,specres=specres,vdisp=vdisp,$
         sigmax=sigmax,zsnrcross=zsnrcross,galaxy=galaxy,outpath=specfitpath,suffix=suffix+'_all',$
         templatefile=templatefile,linefile=linefile,indexfile=indexfile,specfit=specfit_all,$
         /nologfile,vmaxtweak=vmaxtweak,vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,$
         doplot=doplot,specdatafile=specdatafile_all)

; the following objects require two-component balmer-line fits, so
; refit them here: NGC1566, NGC3031, NGC4579, NGC5033; candidate AGN
; that were deemed to *not* have a broad component: NGC1266, NGC1316,
; NGC4125, NGC4450, NGC4594, NGC4736, NGC5195

       broad = 'NGC'+['1566','3031','4579','5033']
       broadindx = speclinefit_locate(sings,broad)
       good = where(broadindx ne -1L,ngood)
       if (ngood ne 0L) then begin
          broadindx = broadindx[good]
          sigmax_broad = replicate(1500.0,ngood)
       endif
; trick MPFIT into fitting the narrow and broad lines appropriately by
; giving an input line-width; the dimensions of SIGGUESS must match
; the number of emission lines in LINEFILE_BROAD!
       sigguess = [replicate(100.0,10),replicate(1000.0,3)] 

       specdata_broad = ispeclinefit(wave[*,broadindx],flux[*,broadindx],invvar[*,broadindx],$
         zobj=zobj[broadindx],specres=specres,vdisp=vdisp[broadindx],sigmax=sigmax_broad,$
         zsnrcross=zsnrcross[broadindx],galaxy=galaxy[broadindx],outpath=specfitpath,$
         suffix=suffix+'_broad',templatefile=templatefile,linefile=linefile_broad,$
         indexfile=indexfile,specfit=specfit_broad,/nologfile,vmaxtweak=vmaxtweak[broadindx],$
         vlinemaxtweak=vlinemaxtweak[broadindx],/clobber,silent=silent,doplot=doplot,$
         specdatafile=specdatafile_broad,sigguess=sigguess)

; now prepend the SINGS structure, merge the 'all' data structure with
; the fits to the objects requiring broad lines, and then write out

       specdata = struct_addtags(specdata_all,replicate(linefit_broad,n_elements(specdata_all)))
       junk = specdata[broadindx] & struct_assign, specdata_broad, junk & specdata[broadindx] = junk

       specfit = specfit_all
       specfit[*,*,broadindx] = specfit_broad
       
       specdatafile = repstr(specdatafile_all,'_all','')
       mwrfits, specdata, specdatafile, /create
       spawn, 'gzip -f '+specdatafile, /sh
       
       specfitfile = repstr(specdatafile,'_specdata','_specfit')
       mwrfits, specfit, specfitfile, /create
       spawn, 'gzip -f '+specfitfile, /sh

       if (not keyword_set(test)) then sings_parse_specfit, /nuclear

    endif
       
; ##################################################    
; fit the 20" scans

    if keyword_set(drift20) then begin

       if keyword_set(test) then suffix = 'test_drift20' else suffix = 'sings_drift20'

       d20 = where(sings1.drift20,nd20)
       sings = sings1[d20]
       if (n_elements(index1) ne 0L) then sings = sings1[d20[index1]]
       nd20 = n_elements(sings)

       speclist = strtrim(sings.drift20_file,2)
       galaxy = strtrim(sings.galaxy,2)

; read the spectra; linearly extrapolate the wavelength vector to
; avoid NaN in ibackfit 

       forage = iforage(spec1dpath+speclist)
       npix = max(forage.naxis1)
       wave = dblarr(npix,nd20)
       flux = wave*0.0D
       invvar = wave*0.0D
       for id20 = 0L, nd20-1L do begin
          spec1 = rd1dspec(speclist[id20],datapath=spec1dpath)
          wave[*,id20] = interpol(spec1.wave,lindgen(spec1.npix),lindgen(npix))
;         wave[0L:spec1.npix-1L,id20] = spec1.wave
          flux[0L:spec1.npix-1L,id20] = spec1.spec
          invvar[0L:spec1.npix-1L,id20] = 1.0/spec1.sigspec^2.0
          icleanup, spec1
       endfor
       
; read in the parameters of the fitting       

       d20params1 = rsex(specfitpath+'zguess_drift20.sex')
       these = speclinefit_locate(sings,d20params1.galaxy)
       good = where(these ne -1,ngood)
       if (ngood ne 0L) then d20params = d20params1[good]

       vdisp = d20params.vdisp
       vmaxtweak = d20params.vmaxtweak
       vlinemaxtweak = d20params.vlinemaxtweak
       sigmax = d20params.sigmax
       zsnrcross = d20params.zsnrcross
       zobj = d20params.zobj
       
       specres = 7.5 ; [FWHM, A]

; now fit them!       

       specdata_all = ispeclinefit(wave,flux,invvar,zobj=zobj,specres=specres,vdisp=vdisp,$
         sigmax=sigmax,zsnrcross=zsnrcross,galaxy=galaxy,outpath=specfitpath,suffix=suffix+'_all',$
         templatefile=templatefile,linefile=linefile,indexfile=indexfile,specfit=specfit_all,$
         /nologfile,vmaxtweak=vmaxtweak,vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,$
         doplot=doplot,specdatafile=specdatafile_all)

; objects requiring broad lines
       broad = 'NGC'+['1566','3031','4579','5033']
       broadindx = speclinefit_locate(sings,broad)
       good = where(broadindx ne -1L,ngood)
       if (ngood ne 0L) then begin
          broadindx = broadindx[good]
          sigmax_broad = replicate(2500.0,ngood) ; the larger sigmax is needed for NGC4579
       endif
; trick MPFIT into fitting the narrow and broad lines appropriately by
; giving an input line-width; the dimensions of SIGGUESS must match
; the number of emission lines in LINEFILE_BROAD!
       sigguess = [replicate(100.0,10),replicate(1000.0,3)] 

       specdata_broad = ispeclinefit(wave[*,broadindx],flux[*,broadindx],invvar[*,broadindx],$
         zobj=zobj[broadindx],specres=specres,vdisp=vdisp[broadindx],sigmax=sigmax_broad,$
         zsnrcross=zsnrcross[broadindx],galaxy=galaxy[broadindx],outpath=specfitpath,$
         suffix=suffix+'_broad',templatefile=templatefile,linefile=linefile_broad,$
         indexfile=indexfile,specfit=specfit_broad,/nologfile,vmaxtweak=vmaxtweak[broadindx],$
         vlinemaxtweak=vlinemaxtweak[broadindx],/clobber,silent=silent,doplot=doplot,$
         specdatafile=specdatafile_broad,sigguess=sigguess)

; now prepend the SINGS structure, merge the 'all' data structure with
; the fits to the objects requiring broad lines, and then write out

       specdata = struct_addtags(specdata_all,replicate(linefit_broad,n_elements(specdata_all)))
       junk = specdata[broadindx] & struct_assign, specdata_broad, junk & specdata[broadindx] = junk

       specfit = specfit_all
       specfit[*,*,broadindx] = specfit_broad
       
       specdatafile = repstr(specdatafile_all,'_all','')
       mwrfits, specdata, specdatafile, /create
       spawn, 'gzip -f '+specdatafile, /sh
       
       specfitfile = repstr(specdatafile,'_specdata','_specfit')
       mwrfits, specfit, specfitfile, /create
       spawn, 'gzip -f '+specfitfile, /sh

       if (not keyword_set(test)) then sings_parse_specfit, /drift20

    endif
       
; ##################################################    
; fit the 56" scans

    if keyword_set(drift56) then begin

       if keyword_set(test) then suffix = 'test_drift56' else suffix = 'sings_drift56'

       d56 = where(sings1.drift56,nd56)
       sings = sings1[d56]
       if (n_elements(index1) ne 0L) then sings = sings1[d56[index1]]
       nd56 = n_elements(sings)

       speclist = strtrim(sings.drift56_file,2)
       galaxy = strtrim(sings.galaxy,2)

; read the spectra; linearly extrapolate the wavelength vector to
; avoid NaN in ibackfit 

       forage = iforage(spec1dpath+speclist)
       npix = max(forage.naxis1)
       wave = dblarr(npix,nd56)
       flux = wave*0.0D
       invvar = wave*0.0D
       for id56 = 0L, nd56-1L do begin
          spec1 = rd1dspec(speclist[id56],datapath=spec1dpath)
          wave[*,id56] = interpol(spec1.wave,lindgen(spec1.npix),lindgen(npix))
;         wave[0L:spec1.npix-1L,id56] = spec1.wave
          flux[0L:spec1.npix-1L,id56] = spec1.spec
          invvar[0L:spec1.npix-1L,id56] = 1.0/spec1.sigspec^2.0
          icleanup, spec1
       endfor
       
; read in the parameters of the fitting       

       d56params1 = rsex(specfitpath+'zguess_drift56.sex')
       these = speclinefit_locate(sings,d56params1.galaxy)
       good = where(these ne -1,ngood)
       if (ngood ne 0L) then d56params = d56params1[good]

       vdisp = d56params.vdisp
       vmaxtweak = d56params.vmaxtweak
       vlinemaxtweak = d56params.vlinemaxtweak
       sigmax = d56params.sigmax
       zsnrcross = d56params.zsnrcross
       zobj = d56params.zobj
       
       specres = 7.5 ; [FWHM, A]

; now fit them!       

       specdata_all = ispeclinefit(wave,flux,invvar,zobj=zobj,specres=specres,vdisp=vdisp,$
         sigmax=sigmax,zsnrcross=zsnrcross,galaxy=galaxy,outpath=specfitpath,suffix=suffix+'_all',$
         templatefile=templatefile,linefile=linefile,indexfile=indexfile,specfit=specfit_all,$
         /nologfile,vmaxtweak=vmaxtweak,vlinemaxtweak=vlinemaxtweak,/clobber,silent=silent,$
         doplot=doplot,specdatafile=specdatafile_all)

; although none of the objects had broad lines, we still need to add
; the relevant structure tags and write out

       specdata = struct_addtags(specdata_all,replicate(linefit_broad,n_elements(specdata_all)))
       specfit = specfit_all
       
       specdatafile = repstr(specdatafile_all,'_all','')
       mwrfits, specdata, specdatafile, /create
       spawn, 'gzip -f '+specdatafile, /sh
       
       specfitfile = repstr(specdatafile,'_specdata','_specfit')
       mwrfits, specfit, specfitfile, /create
       spawn, 'gzip -f '+specfitfile, /sh

       if (not keyword_set(test)) then sings_parse_specfit, /drift56

    endif
       
return
end    

;;; ##################################################    
;;; fit the 20" scans; increase VMAXSHIFT to 400 km/s; I have
;;; verified (jm07aug07nyu) that this does not adversely affect the
;;; continuum fits on objects with strong Balmer absorption
;;
;;    if keyword_set(drift20) then begin
;;
;;       if keyword_set(test) then suffix = 'test' else suffix = 'sings_drift20'
;;
;;       if (n_elements(galaxy) eq 0L) then $
;;         d20 = where(sings.drift20 and (sings.drift20_agnflag eq 0L)) else $
;;           d20 = speclinefit_locate(sings,galaxy)
;;       speclist = strtrim(sings[d20].drift20_file,2)
;;       nd20 = n_elements(d20)
;;
;;       specres = fltarr(nd20)+7.9 ; KPNO
;;       ctio = where(strmatch(sings[d20].drift20_observat,'*ctio*',/fold))
;;       if (ctio[0] ne -1L) then specres[ctio] = 8.4
;;
;;;      linefile = 'elinelist_drift20.dat'
;;
;;       specdata = sings_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,$
;;         datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
;;         /charlot,/zcrosscor,/postscript,/write,vmaxshift=500.0,/nologfile,$
;;         starvdisp=100.0,/drift20,/positive_emission,specdatafile=specdatafile,$
;;         eigendir=eigendir,eigenfile=eigenfile,linefile=linefile)
;;
;;       if (not keyword_set(test)) then sings_parse_specfit, /drift20
;;
;;    endif
;;
;;; ##################################################    
;;; fit the 56" scans; increase VMAXSHIFT to 350 km/s; I have verified
;;; (jm07aug07nyu) that this does not adversely affect the continuum
;;; fits on objects with strong Balmer absorption
;;
;;    if keyword_set(drift56) then begin
;;
;;       if keyword_set(test) then suffix = 'test' else suffix = 'sings_drift56'
;;
;;       if (n_elements(galaxy) eq 0L) then $
;;         d56 = where(sings.drift56 and (sings.drift56_agnflag eq 0L)) else $
;;           d56 = speclinefit_locate(sings,galaxy)
;;       speclist = strtrim(sings[d56].drift56_file,2)
;;       nd56 = n_elements(d56)
;;
;;       specres = fltarr(nd56)+7.9 ; KPNO
;;       ctio = where(strmatch(sings[d56].drift56_observat,'*ctio*',/fold))
;;       if (ctio[0] ne -1L) then specres[ctio] = 8.4
;;       
;;;      linefile = 'elinelist_drift56.dat'
;;
;;       specdata = sings_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,$
;;         datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
;;         /charlot,/zcrosscor,/postscript,/write,vmaxshift=500.0,/nologfile,$
;;         starvdisp=100.0,/drift56,/positive_emission,specdatafile=specdatafile,$
;;         eigendir=eigendir,eigenfile=eigenfile,linefile=linefile)
;;
;;       if (not keyword_set(test)) then sings_parse_specfit, /drift56
;;
;;    endif
;;
;;stop    
    
; replace some of the measured fluxes and EWs with upper limits: these
; objects and lines were individually examined (jm07aug08nyu) and were
; usually related to poor continuum subtraction and problems with
; [SII] in the noisy red part of the spectra

;      old = mrdfits(specfitpath+specdatafile+'.gz',1,/silent)
;      new = old
;
;      indx = speclinefit_locate(new,'ngc1404')
;      new[indx].h_gamma     = new[indx].h_gamma_limit
;      new[indx].h_beta      = new[indx].h_beta_limit
;      new[indx].h_alpha     = new[indx].h_alpha_limit
;      new[indx].sii_6716    = new[indx].sii_6716_limit
;      new[indx].h_gamma_ew  = new[indx].h_gamma_ew_limit
;      new[indx].h_beta_ew   = new[indx].h_beta_ew_limit
;      new[indx].h_alpha_ew  = new[indx].h_alpha_ew_limit
;      new[indx].sii_6716_ew = new[indx].sii_6716_ew_limit
;      
; old code:
;      specdata = sings_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,$
;        datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
;        /charlot,/zcrosscor,/postscript,/write,vmaxshift=300.0,/nologfile,$
;        starvdisp=100.0,/drift20,/positive_emission,specdatafile=specdatafile,$
;        eigendir=eigendir,eigenfile=eigenfile,linefile=linefile)

; old code:
;      specdata = sings_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,$
;        datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
;        /charlot,/zcrosscor,/postscript,/write,vmaxshift=300.0,/nologfile,$
;        starvdisp=100.0,/nuclear,/positive_emission,specdatafile=specdatafile,$
;        eigendir=eigendir,eigenfile=eigenfile,linefile=linefile)

; old code:
;      specdata = sings_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,$
;        datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
;        /charlot,/zcrosscor,/postscript,/write,vmaxshift=300.0,/nologfile,$
;        starvdisp=100.0,/drift56,/positive_emission,specdatafile=specdatafile,$
;        eigendir=eigendir,eigenfile=eigenfile,linefile=linefile)

