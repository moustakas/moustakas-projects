;+
; NAME:
;       ISPECLINEFIT()
;
; PURPOSE:
;       Fit a list of galaxy spectra.
;
; INPUTS:
;       speclist - FITS file list of galaxy spectra to fit [NSPEC]
;       specres  - instrumental resolution [FWHM, Angstrom]; SPECRES
;                  can be a scalar or a [2,NSPECRES] array, where the
;                  zeroth dimension is the wavelength and the first
;                  dimension contains the spectral resolution at each
;                  wavelength 
;
; OPTIONAL INPUTS:
;       datapath  - I/O path
;       eigenspec - pre-defined spectral template sets (also see
;                   CHABRIER, POLYNOMIALS, and EIGENFILE)
;          0 - BC03, Z = 0.02 (default)
;          1 - BC03, Z = 0.004
;          2 - BC03, Z = 0.05
;       eigenfile - fitting template set (overwrites EIGENSPEC) 
;       eigendir  - path to the pre-defined templates in EIGENSPEC or
;                   to EIGENFILE (default ${ISPEC_DIR}/templates)
;       linefile  - emission line file (default 'elinelist.dat')  
;       linepath  - path to LINEFILE (default DATAPATH)
;       dustmodel - continuum dust model
;          0 - tied reddening for all templates (default)
;          1 - unique reddening for all templates
;          2 - templates older than 30 Myr have reddening that is
;              [0.1,0.9] times the reddening of the young templates 
;          3 - old-to-young reddening ratio ("alpha") is fixed at 0.3 
;       snrcut    - compute upper limits on lines with S/N < SNRCUT
;                   (default 1.0)
;       npoly     - number of polynomial templates if POLYNOMIALS=1
;                   (default 10) 
;       psname    - postscript output file name if POSTSCRIPT=1
;                   (default MJDSTR+'_'+SUFFIX+'_specfit.ps', where MJDSTR
;                   is the modified Julian date)
;       suffix    - unique output file name suffix (default '')
;       starvdisp - overwrite the velocity dispersion VDISP in the
;                   header
;       zobj      - galaxy redshift [NSPEC]
;       extra     - keywords and parameters for IFITSPEC
;
; KEYWORD PARAMETERS:
;       Zmulti      - find the best-fitting stellar metallicity
;                     templates from a pre-defined set of three
;                     template sets (stronger than EIGENSPEC and
;                     EIGENFILE) 
;       chabrier    - use the Chabrier IMF (default Salpeter IMF) 
;       polynomials - fit a linear combination of NPOLY polynomials
;                     (stronger than EIGENSPEC, EIGENFILE, and
;                     ZMULTI); also sets the following keywords:
;                     NBACK=0, NDUST=0, and BACKNEGATIVE=1)  
;       zcrosscor   - use cross-correlation to update the
;                     absorption-line redshift (see IFITSPEC) 
;       doplot      - generate screen plots (if POSTSCRIPT=1 then
;                     DOPLOT=0)   
;       debug       - enable additional diagnostic plots in IFITSPEC
;                     (if POSTSCRIPT=1 then DEBUG=0) 
;       postscript  - generate postscript output of the fitting
;                     results (see PSNAME)
;       write       - write the results of the spectral line and
;                     continuum fitting as two binary FITS files (see
;                     COMMENTS)  
;       nologfile   - do not generate a log file
;
; OUTPUTS:
;       specdata   - line-fitting results (fluxes, EWs, widths, etc.)
;
; OPTIONAL OUTPUTS:
;       mjdstr     - modified Julian date prefix string
;       specfit    - fitting results for the last object, useful for
;                    debugging ([4,NPIX] array)
;          specfit[0,*] = rest wavelength
;          specfit[1,*] = rest flux
;          specfit[2,*] = best-fitting continuum spectrum
;          specfit[3,*] = best-fitting emission-line spectrum 
;
; COMMENTS:
;
; PROCEDURES USED:
;       CWD(), SPLOG, MRDFITS(), MAKE_WAVE(), SXPAR(), MWRFITS,
;       IM_OPENCLOSE, READ_LINEPARS(), GET_JULDATE, RD1DSPEC(),
;       STRUCT_TRIMTAGS(), IFITSPEC, STRUCT_ADDTAGS(),
;       PARSE_ILINEFIT(), MKHDR, SXDELPAR, SXADDPAR,
;       IFIGURE_SPECLINEFIT, ICLEANUP, LINTERP, VARIANCE()
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003-2004, U of A - written
;       jm04mar03uofa - added ZMULTI, CHABRIER, and DUSTMODEL keywords 
;       jm04apr29uofa - added CONTINUUM_POLYCOEFF structure tag
;       jm04may09uofa - added STARVDISP input to IFITSPEC and to
;                       SPECDATA output structure; add the STEP field
;                       to the ALPHA structure
;       jm04jun16uofa - overwrite the default STARVDISP if requested 
;       jm04jul23uofa - lower limit velocity dispersion = 50 km/s
;       jm05jan26uofa - added NOLOGFILE keyword
;       jm05jul22uofa - removed ZOBJ_ERR from the output data
;                       structure; no longer read from the header 
;       jm05jul26uofa - added ZOBJ optional input
;       jm07sep23nyu - better error handling of SPECRES; also SPECRES
;                      can now be different for each object
;
; Copyright (C) 2003-2005, 2007, John Moustakas
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

function atlas_ispeclinefit, speclist, specres=specres1, datapath=datapath, eigenspec=eigenspec, $
  eigenfile=eigenfile, eigendir=eigendir, linefile=linefile, linepath=linepath, $
  dustmodel=dustmodel, snrcut=snrcut, npoly=npoly, psname=psname, suffix=suffix, $
  starvdisp=starvdisp1, zobj=zobj1, _extra=extra, mjdstr=mjdstr, specfit=specfit, $
  Zmulti=Zmulti, chabrier=chabrier, polynomials=polynomials, zcrosscor=zcrosscor, $
  doplot=doplot, debug=debug, postscript=postscript, write=write, nologfile=nologfile, $
  integrated=integrated, nuclear=nuclear

    nspec = n_elements(speclist)
    if (nspec eq 0L) then begin
       doc_library, 'ispeclinefit'
       return, -1 
    endif

    light = 2.99792458D5 ; speed of light [km/s]

    if n_elements(datapath) eq 0L then datapath = cwd()
    if n_elements(snrcut) eq 0L then snrcut = 1.0
    if n_elements(suffix) eq 0L then suffix = ''

    nzobj = n_elements(zobj1)
    if (nzobj ne 0L) then if (nzobj ne nspec) then begin
       splog, 'ZOBJ must be the same size as SPECLIST!'
       return, -1L
    endif

; pre-defined templates

    if n_elements(eigendir) eq 0L then eigendir = filepath('',$
      root_dir=getenv('ISPEC_DIR'),subdirectory='templates')

    if keyword_set(chabrier) then imf = 'chabrier' else imf = 'salpeter'

    if n_elements(eigenspec) eq 0L then eigenspec = 0L
    if (n_elements(eigenfile) eq 0L) then begin

       case eigenspec of
          1L:   eigenfile = 'BC03_Z004_'+imf+'_templates.fits' ; LMC metallicity 
          2L:   eigenfile = 'BC03_Z05_'+imf+'_templates.fits'  ; twice-solar metallicity 
          else: eigenfile = 'BC03_Z02_'+imf+'_templates.fits'  ; solar metallicity 
       endcase 

    endif 

; pure polynomial template set; particular IFITSPEC keywords must be
; overwritten 
    
    if keyword_set(polynomials) then begin

       if n_elements(npoly) eq 0L then npoly = 10L

       eigenfile = 'Polynomials ('+string(npoly,format='(I0)')+')'
       eigeninfo = {$
         ntemplate:             npoly, $
         template_name: 'Polynomials', $
         template_age:  fltarr(npoly), $
         template_Z:              0.0, $   
         template_imf:          'N/A', $
         template_ML_U: fltarr(npoly)+1.0, $
         template_ML_B: fltarr(npoly)+1.0, $
         template_ML_V: fltarr(npoly)+1.0, $
         template_ML_R: fltarr(npoly)+1.0}
;        template_ML_I: fltarr(npoly)+1.0}
       
       eigen_minwave = 2500.0
       eigen_maxwave = 25000.0
       eigen_cd1_1 = 5.0
       
       eigen_npix = fix((eigen_maxwave-eigen_minwave)/eigen_cd1_1)+1
       eigenwave = findgen(eigen_npix)*eigen_cd1_1+eigen_minwave
       eigenflux = poly_array(eigen_npix,npoly)
       eigenres = -1.0

       zcrosscor = 0L
       dustmodel = 0L

       if (n_elements(extra) eq 0L) then begin
          extra = {nback: 0L, ndust: 0L, backnegative: 1L}
       endif else begin
          if tag_exist(extra,'NBACK') then extra.nback = 0L else extra = create_struct(extra,'nback',0L)
          if tag_exist(extra,'NDUST') then extra.ndust = 0L else extra = create_struct(extra,'ndust',0L)
          if tag_exist(extra,'BACKNEGATIVE') then extra.backnegative = 1L else $
            extra = create_struct(extra,'backnegative',1L)
       endelse 
             
    endif

; read the eigen-templates

    if (not keyword_set(polynomials)) then begin
    
; multi-metallicity template set (pre-defined)
       
       if keyword_set(Zmulti) then begin

          splog, 'ZMULTI keyword set.'
          bigeigenfile = [$
            'BC03_Z004_'+imf+'_templates.fits', $
            'BC03_Z02_'+imf+'_templates.fits', $
            'BC03_Z05_'+imf+'_templates.fits']

          for i = 0L, n_elements(bigeigenfile)-1L do begin

             if file_test(eigendir+bigeigenfile[i],/regular) eq 0L then begin
                splog, 'Eigen-template file '+eigendir+bigeigenfile[i]+' not found.'
                return, -1L
             endif

             splog, 'Reading '+eigendir+bigeigenfile[i]+'.'
             eigeninfo1 = mrdfits(eigendir+bigeigenfile[i],1,/silent)           ; info structure
             eigenflux1 = mrdfits(eigendir+bigeigenfile[i],2,eigenhead,/silent) ; flux

             if (i eq 0L) then begin
             
                eigenwave = make_wave(eigenhead)       ; wavelength
                eigenres = sxpar(eigenhead,'EIGENRES') ; FWHM resolution [Angstrom] 

                eigeninfo = eigeninfo1
                eigenflux = eigenflux1
                
             endif else begin

                eigeninfo = [ [eigeninfo], [eigeninfo1] ]
                eigenflux = [ [ [eigenflux] ], [ [eigenflux1] ] ]
                
             endelse
             
          endfor

          eigeninfo = reform(eigeninfo)
          eigenflux = reform(eigenflux)

       endif else begin 
          
          if file_test(eigendir+eigenfile,/regular) eq 0L then begin
             splog, 'Eigen-template file '+eigendir+eigenfile+' not found.'
             return, -1L
          endif

          splog, 'Reading '+eigendir+eigenfile+'.'
          eigeninfo = mrdfits(eigendir+eigenfile,1,/silent)           ; info structure
          eigenflux = mrdfits(eigendir+eigenfile,2,eigenhead,/silent) ; flux
          eigenwave = make_wave(eigenhead)                            ; wavelength
          eigenres = sxpar(eigenhead,'EIGENRES')                      ; FWHM resolution [Angstrom] 
          
       endelse
          
; consider the DUSTMODEL keyword

       if n_elements(dustmodel) eq 0L then dustmodel = 0L

       splog, 'DUSTMODEL = '+string(dustmodel,format='(I0)')
       alpha = {$
         indx:             0L, $
         fixed:            1L, $
         value:          1.0D, $
         limited:     [0L,0L], $
         limits:  [1.0D,1.0D]  $
         }
       alpha = replicate(alpha,eigeninfo[0].ntemplate)

       case dustmodel of
          0L: begin ; alpha and reddening fixed to the same value
             redindx = lonarr(eigeninfo[0].ntemplate)
             alpha.indx = lonarr(eigeninfo[0].ntemplate)
          end
          1L: begin ; alpha and reddening free for every template
             redindx = lindgen(eigeninfo[0].ntemplate)
             alpha.indx = lindgen(eigeninfo[0].ntemplate)
          end
          2L: begin ; alpha=[0.1,0.9] for old stellar pops
             young = where(eigeninfo[0].template_age le 1D7,nyoung,comp=old,ncomp=nold)
             redindx = lonarr(eigeninfo[0].ntemplate)
             alpha[old].indx = 1L    ; tie all the old templates to one another
             alpha[old].fixed = 0L
             alpha[old].value = 0.3D ; starting guess
             alpha[old].limited = 1L
;            alpha[old].limits = [0.1D,0.9D]
             alpha[old].limits = [0.0D,1.0D]
          end
          3L: begin ; alpha=0.3 for old stellar pops
             young = where(eigeninfo[0].template_age le 1D7,nyoung,comp=old,ncomp=nold)
             redindx = lonarr(eigeninfo[0].ntemplate)
             alpha[old].indx = 1L ; tie all the old templates to one another
             alpha[old].value = 0.3D
          end
          else: redindx = lonarr(eigeninfo[0].ntemplate)
       endcase
             
    endif 

    best_eigeninfo = eigeninfo ; we need this for the loop below

; read the emission line parameters and constraints

    if n_elements(linefile) eq 0L then linefile = 'elinelist.dat'
    if n_elements(linepath) eq 0L then linepath = datapath

    linepars = read_linepars(linepath=linepath,linefile=linefile)
    if size(linepars,/type) ne 8L then return, -1L
    nline = n_elements(linepars)

    zindex = linepars.zindex
    windex = linepars.windex

; error checking on the spectral resolution vector; if SPECRES is a
; scalar then the square brackets ensures that SPECRES_NDIM and
; SPECRES_DIM are unity

    specres_dim = size([specres1],/dimension)
    specres_ndim = size([specres1],/n_dimension)

    case specres_ndim of
       0L: begin
          splog, 'SPECRES must be specified.'
          return, -1L
       end
       1L: begin
          nspecres = size(specres1,/n_elements)      
          if (nspecres eq 1L) then specres = replicate(specres1,nspec) else begin
             if (nspecres ne nspec) then begin
                splog, 'Dimensions of SPECRES and SPECLIST must agree.'
                return, -1L
             endif else specres = specres1
          endelse
       end
       2L: ; handle below
       else: begin
          splog, 'SPECRES must be a one- or two-dimensional array.'
          return, -1L
       end
    endcase

; open the log file

    get_juldate, jd
    mjdstr = string(long(jd-2400000L),format='(I5)')

    if keyword_set(write) then begin

       if (suffix ne '') then logfile = mjdstr+'_'+suffix+'.log' else logfile = mjdstr+'.log'
       if (suffix ne '') then $
         specfitfile = mjdstr+'_'+suffix+'_specfit.fits' else $
         specfitfile = mjdstr+'_specfit.fits'
       if (suffix ne '') then $
         specdatafile = mjdstr+'_'+suffix+'_specdata.fits' else $
         specdatafile = mjdstr+'_specdata.fits'

       if (not keyword_set(nologfile)) then begin
          splog, filename=logfile
          splog, 'Log file '+logfile+' opened '+systime()
       endif

    endif

    splog, 'IDL version: ' + string(!version,format='(99(A," "))')
    spawn, ['uname -a'], uname
    splog, 'UNAME: '+uname[0]
    
    splog, 'Current datapath is ', datapath
    splog, 'Fitting with '+string(eigeninfo[0].ntemplate,format='(I0)')+' templates.'
    
; we will store the results of the continuum fitting for each galaxy
; in SPECLIST in a single extension of ?????_specfit.fits.
; consequently in the zeroth extension we generate a table of contents
; (toc) structure describing the contents of each extension

; table of contents structure for ?????_SPECFIT.FITS

    toc = {$
      fitdate:   im_today(), $
      mjdstr:    mjdstr, $
      nfitspec:  nspec, $
      ap00:      'rest wavelength', $
      ap01:      'rest data spectrum', $
      ap02:      'rest continuum spectrum', $
      ap03:      'rest emission line spectrum'}
    
    if keyword_set(write) then begin

       splog, 'Creating '+specfitfile+'.'
       mwrfits, toc, specfitfile, /create

    endif

    if n_elements(psname) eq 0L then if (suffix ne '') then $
      psname = mjdstr+'_'+suffix+'_specfit' else $
      psname = mjdstr+'_specfit'

    if keyword_set(postscript) then begin
       doplot = 0
       debug = 0
       splog, 'Opening plot file '+psname
       im_openclose, psname, /courier, /postscript, /silent, /landscape
    endif
 
; ---------------------------------------------------------------------------
; loop on each spectrum in SPECLIST
; ---------------------------------------------------------------------------
    
    stime0 = systime(1)
    for iobj = 0L, nspec-1L do begin

       splog, format='("Fitting object #",I0,"/",I0,".")', iobj+1, nspec
       
; read the spectrum and the redshift

       scube = rd1dspec(speclist[iobj],datapath=datapath)
       flux = scube.spec
       ferr = scube.sigspec
       wave = scube.wave
       header = scube.header
       npix = scube.npix

       invvar = flux*0D0+1.0/(djs_median(flux)/50.0)^2.0 ; "intelligent" default INVVAR
       good = where(ferr gt 0.0,ngood)
       if (ngood gt 0L) then begin
          invvar[good] = 1.0/double(ferr[good])^2.0
          median_snr = djs_median(flux[good]/ferr[good])
       endif else begin
          splog, 'WARNING: Improper error spectrum -- errors will be wrong!'
          median_snr = -1.0
       endelse
       
; instrumental resolution; resample, do not interpolate between pixels  

       case specres_ndim of
          1L: fitspecres = specres[iobj]
          2L: linterp, reform(specres[0,*]), reform(specres[1,*]), wave, fitspecres, /nointerp
       endcase

; extract redshift information on GALAXY
       
       galaxy = strcompress(sxpar(header,'GALAXY'),/remove)
       if (nzobj ne 0L) then zobj = zobj1[iobj] else begin
          zobj = float(sxpar(header,'Z',count=zcount))
          if (zcount eq 0L) then begin
             splog, 'WARNING: No redshift information for '+strtrim(speclist[iobj],2)+' (assuming z=0.0).'
             zobj = 0.0
          endif
       endelse

; ###########################################################################
; SPECIAL CASES; jm06dec07nyu

       no_medcontinuum = 0L ; default is to median-subtract the residuals
       medsmooth_window = 150L  ; larger than the default values!
       smooth_window = 50L
       zcrosscor = 1L
       
       linepars = read_linepars(linepath=linepath,linefile=linefile)
       zindex = linepars.zindex & windex = linepars.windex
       zline = zobj
       sigmax = 500.0

; INTEGRATED SPECTRA

       if keyword_set(integrated) then begin

          if strmatch(galaxy,'*mcg-03-04-014*',/fold) then zobj = 0.036
          if strmatch(galaxy,'*ic0860*',/fold) then zobj = 0.0135
          if strmatch(galaxy,'*ngc1800*',/fold) then zobj = 0.0042
          if strmatch(galaxy,'*ngc7713*',/fold) then zobj = 0.0039
          if strmatch(galaxy,'*ngc0232*',/fold) then zobj = 0.024
          if strmatch(galaxy,'*ngc5264*',/fold) then zobj = 0.004
          if strmatch(galaxy,'*ngc7130*',/fold) then zobj = 0.018
          if strmatch(galaxy,'*ic5179*',/fold) then zobj = 0.013
          if strmatch(galaxy,'*ngc3991n',/fold) then zcrosscor = 0L ; do not cross-correlate!
          if strmatch(galaxy,'*ugca116*',/fold) then no_medcontinuum = 1L

; separately fit [SII] and [NII] line-width/redshift
          if $
            strmatch(galaxy,'*ngc0922*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_integrated_ngc0922.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; separately fit the [SII] line-width/redshift
          if $
            strmatch(galaxy,'*ngc1275*',/fold) or $
            strmatch(galaxy,'*ngc4676*',/fold) or $ ; = NGC4676, NGC4676A, and NGC4676B
            strmatch(galaxy,'*ic5298*',/fold) or $
            strmatch(galaxy,'*ngc5256*',/fold) or $
            strmatch(galaxy,'*ngc7713*',/fold) or $
            strmatch(galaxy,'*ngc7771*',/fold) or $
            strmatch(galaxy,'*ugc12747*',/fold) or $
            strmatch(galaxy,'*ngc7591*',/fold) or $
            strmatch(galaxy,'*ngc7592*',/fold) or $ ; = NGC7592, NGC7592A, and NGC7592B
            strmatch(galaxy,'*ngc6240*',/fold) or $
            strmatch(galaxy,'*ngc7130*',/fold) or $
            strmatch(galaxy,'*ngc7674*',/fold) or $
            strmatch(galaxy,'*arp182*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_sii_separate.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; the [NII] doublet ratio has to be relaxed, and the [OIII] doublet
; has to be fitted separately
          if $
            strmatch(galaxy,'*iras05189-2524*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_integrated_iras05189-2524.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; broad H-alpha, but not the other Balmer lines
          if $
            strmatch(galaxy,'*ngc7469*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_integrated_ngc7469.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; broad H-alpha, but not the other Balmer lines
          if $
            strmatch(galaxy,'*mrk0315*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_integrated_mrk0315.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; custom care: [NII] and [SII] have to be separated, and I need a
; broad H-alpha line
          if $
            strmatch(galaxy,'*ngc1068*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_integrated_ngc1068.dat')
             zindex = linepars.zindex & windex = linepars.windex
             sigmax = 1500.0
          endif

       endif

; NUCLEAR SPECTRA       
       
       if keyword_set(nuclear) then begin

; NGC1068 needs some hand-tweaking and the fit still isn't
; great, but good enough for goverment work (jm07sep27nyu)           
          if $
            strmatch(galaxy,'*ngc1068*',/fold) then begin
             medsmooth_window = 350L ; larger than the default values!
             smooth_window = 200L
             maskwidth = 100L
             zcrosscor = 0L
             zobj = 0.003597
          endif

; test code for NGC1275 (jm07sep27nyu); needs more work
;         if $
;           strmatch(galaxy,'*ngc1275*',/fold) then begin
;            medsmooth_window = 350L ; larger than the default values!
;            smooth_window = 200L
;            maskwidth = 100L
;            zcrosscor = 0L
;            zobj = 0.01682
;         endif

; separately fit the [SII] line-width
          if $
            strmatch(galaxy,'*ugca166*',/fold) or $
            strmatch(galaxy,'*ngc2623*',/fold) or $
            strmatch(galaxy,'*ngc4051*',/fold) or $
            strmatch(galaxy,'*ugc04459*',/fold) or $
            strmatch(galaxy,'*ngc3998*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_sii_separate.dat')
             zindex = linepars.zindex & windex = linepars.windex
          endif

; the following objects require two-component balmer-line fits

          if $
            strmatch(galaxy,'*mrk0315*',/fold) then begin
             linepars = read_linepars(linepath=linepath,linefile='elinelist_nuclear_mrk0315.dat')
             zindex = linepars.zindex & windex = linepars.windex
             sigmax = 1500.0
          endif

;; the following objects require two-component balmer-line fits
;
;          if $
;            strmatch(galaxy,'*mrk0315*',/fold) then begin
;             linepars = read_linepars(linepath=linepath,linefile='elinelist_balmer_broad.dat')
;             zindex = linepars.zindex & windex = linepars.windex
;             sigmax = 1500.0
;          endif

       endif
; ###########################################################################
          
       splog, 'Redshift = '+strtrim(string(zobj,format='(F11.6)'),2)+'.'

; extract rest-frame velocity dispersion information on GALAXY

       if (n_elements(starvdisp1) eq 0L) then begin

          starvdisp = float(sxpar(header,'VDISP',count=vcount))
          if (vcount eq 0L) then begin
             starvdisp = 100.0  ; [km/s]
             splog, 'WARNING: No stellar velocity dispersion for '+$
               strtrim(speclist[iobj],2)+'.'
          endif else if (starvdisp le 50.0) then begin
             starvdisp = 50.0  ; [km/s]
             splog, 'WARNING: Stellar velocity dispersion <50 km/s for '+$
               strtrim(speclist[iobj],2)+'.'
          endif

       endif else starvdisp = starvdisp1
       splog, 'Velocity dispersion = '+strtrim(string(starvdisp,format='(F10.1)'),2)+' km/s.'
          
; fit the full spectrum

       ifitspec, flux, wave, eigenflux, eigenwave, eigenres=eigenres, specres=fitspecres, $
         linepars=linepars, invvar=invvar, starvdisp=starvdisp, zobj=zobj, zline=zline, snrcut=snrcut, $
         zindex=zindex, windex=windex, findex=findex, fvalue=fvalue, niter=2L, $ ; note!!$
         _extra=extra, backfit=backfit, linefit=linefit, balmerabs=balmerabs, indices=indices, $
         speclinefit=speclinefit, starflux=starflux, Zbestindx=Zbestindx, debug=debug, $
         zcrosscor=zcrosscor, /combine_blends, Zmulti=Zmulti, redindx=redindx, alpha=alpha, $
         medsmooth_window=medsmooth_window, smooth_window=smooth_window, maskwidth=maskwidth, $
         no_medcontinuum=no_medcontinuum

       if keyword_set(Zmulti) then best_eigeninfo = eigeninfo[Zbestindx]
       
; add broad-line components (customized to SINGS_ISPECLINEFIT!); the
; order of the broad emission lines must match
; ELINELIST_BROAD_NARROW.DAT; the broad emission lines must be sorted
; by decreasing wavelength!!

       if (total(strmatch(linefit.linename,'*h_alpha_broad*',/fold)) eq 0.0) then begin
          linefit_broad = icreate_linefit(1L)
          linefit_broad.linename = 'H_alpha_broad'
          linefit_broad.linewave = 6562.80
          linefit = [linefit,linefit_broad]
       endif
       if (total(strmatch(linefit.linename,'*h_beta_broad*',/fold)) eq 0.0) then begin
          linefit_broad = icreate_linefit(1L)
          linefit_broad.linename = 'H_beta_broad'
          linefit_broad.linewave = 4861.325
          linefit = [linefit,linefit_broad]
       endif
       if (total(strmatch(linefit.linename,'*h_gamma_broad*',/fold)) eq 0.0) then begin
          linefit_broad = icreate_linefit(1L)
          linefit_broad.linename = 'H_gamma_broad'
          linefit_broad.linewave = 4340.464
          linefit = [linefit,linefit_broad]
       endif
       if (total(strmatch(linefit.linename,'*h_delta_broad*',/fold)) eq 0.0) then begin
          linefit_broad = icreate_linefit(1L)
          linefit_broad.linename = 'H_delta_broad'
          linefit_broad.linewave = 4101.734
          linefit = [linefit,linefit_broad]
       endif

       broad = where(strmatch(linefit.linename,'*broad*',/fold),nbroad)
       if (nbroad ne 0L) then begin
          wsort = reverse(sort(linefit[broad].linewave))
          linefit[broad] = linefit[broad[wsort]]
       endif

; parse the results: append the basic data, the absorption-line
; measurements, Lick indices, and continuum-fitting results and
; emission-line measurements together

       outinfo = struct_addtags(struct_addtags(struct_trimtags(toc,except=['AP*']),$
         struct_trimtags(best_eigeninfo,except='TEMPLATE_Z')),$
         {specfile: speclist[iobj], galaxy: galaxy, fit_id: iobj, $
         continuum_minwave: min(wave), continuum_maxwave: max(wave), $
         continuum_snr: float(backfit.continuum_snr), continuum_chi2: float(backfit.continuum_chi2), $
         continuum_vdisp: float(starvdisp), continuum_ebv: float(backfit.starebv), $
         continuum_coeff: double(backfit.starcoeff), continuum_alpha: float(backfit.staralpha), $
         continuum_Z: float(best_eigeninfo.template_Z), continuum_polycoeff: double(backfit.polycoeff), $
         mpfit_status: backfit.mpfit_status, mpfit_niter: backfit.mpfit_niter, $
         z_obj: zobj, z_abs: backfit.z_abs, z_abs_err: backfit.z_abs_err, $
         linefit_chi2: backfit.linefit_chi2, linefit_dof: backfit.linefit_dof, $
         linefit_niter: backfit.linefit_niter, linefit_status: backfit.linefit_status})

       if (n_elements(balmerabs) ne 0L) then begin
          
          babs = create_struct($
            'BABS_'+balmerabs[0].babs_line+'_WAVE',balmerabs[0].babs_wave,$
            'BABS_'+balmerabs[0].babs_line+'_LLINE',balmerabs[0].babs_lline,$
            'BABS_'+balmerabs[0].babs_line+'_ULINE',balmerabs[0].babs_uline,$
            'BABS_'+balmerabs[0].babs_line+'_CONTINUUM',[balmerabs[0].babs_continuum,balmerabs[0].babs_continuum_err],$
            'BABS_'+balmerabs[0].babs_line,[balmerabs[0].babs,balmerabs[0].babs_err],$
            'BABS_'+balmerabs[0].babs_line+'_EW',[balmerabs[0].babs_ew,balmerabs[0].babs_ew_err])
          for k = 1L, n_elements(balmerabs)-1L do babs = create_struct(babs,$
            'BABS_'+balmerabs[k].babs_line+'_WAVE',balmerabs[k].babs_wave,$
            'BABS_'+balmerabs[k].babs_line+'_LLINE',balmerabs[k].babs_lline,$
            'BABS_'+balmerabs[k].babs_line+'_ULINE',balmerabs[k].babs_uline,$
            'BABS_'+balmerabs[k].babs_line+'_CONTINUUM',[balmerabs[k].babs_continuum,balmerabs[k].babs_continuum_err],$
            'BABS_'+balmerabs[k].babs_line,[balmerabs[k].babs,balmerabs[k].babs_err],$
            'BABS_'+balmerabs[k].babs_line+'_EW',[balmerabs[k].babs_ew,balmerabs[k].babs_ew_err])

          babs = struct_addtags({babs_linename: balmerabs.babs_line},babs)
          moreinfo = struct_addtags(babs,indices)

       endif else moreinfo = indices

; parse the LINEFIT structure and append everything together

       plinefit = parse_ilinefit(linefit)
       specdata1 = struct_addtags(struct_addtags(outinfo,{z_line: plinefit.z_line,z_line_err: $
         plinefit.z_line_err}),struct_addtags(moreinfo,struct_trimtags(plinefit,except='Z_LINE*')))       
          
; initialize the output FITS vector and header; everything here is in
; *rest* wavelengths!!

       zplus1 = (specdata1.z_abs+1.0)
       specfit = float([ [wave/zplus1], [zplus1*flux], [backfit.continuum], [speclinefit], [backfit.medcontinuum] ])
;      specfit = float([ [wave/zplus1], [zplus1*flux], [backfit.continuum], [speclinefit] ])

; generate a QA plot of the results

       if keyword_set(doplot) or keyword_set(postscript) then begin

          if keyword_set(postscript) then splog, 'Generating postscript for '+strtrim(speclist[iobj],2)+'.'

          ifigure_speclinefit, specdata1, specfit, starflux, _extra=extra, postscript=postscript

          if (keyword_set(doplot) and (nspec gt 1L)) then begin
             splog, 'Press any key to continue.'
             cc = get_kbrd(1)
          endif

       endif

       if (n_elements(specdata) eq 0L) then specdata = specdata1 else specdata = [ [specdata], [specdata1] ]
       
       if keyword_set(write) then begin

          splog, 'Updating '+specfitfile+'.'
          mwrfits, specfit, specfitfile

; continuously overwrite SPECDATAFILE in case the code crashes
          
          splog, 'Updating '+specdatafile+'.'
          mwrfits, reform(specdata), specdatafile, /create 

       endif 

; clean up memory

       babs = 0 & balmerabs = 0 & backfit = 0 & backfit = 0 & backinfo = 0 
       basicinfo = 0 & flux = 0 & ferr = 0 & invvar = 0 & linefit = 0 
       plinefit = 0 & specdata1 = 0 & indices = 0 & speclinefit = 0
       fitspecres = 0 & wave = 0 & specfit = 0

       icleanup, scube

    endfor 

    splog, format='("Total time for ISPECLINEFIT = ",G0," minutes.")', (systime(1)-stime0)/60.0
    splog, /close
    
; ---------------------------------------------------------------------------
; end main loop
; ---------------------------------------------------------------------------

; close and GZIP files
    
    if keyword_set(write) then begin
       spawn, ['gzip -f '+specdatafile], /sh
       spawn, ['gzip -f '+specfitfile], /sh
    endif

    if keyword_set(postscript) then begin
       im_openclose, postscript=postscript, /close
       spawn, ['gzip -f '+psname+'.ps'], /sh
    endif

return, reform(specdata)
end
