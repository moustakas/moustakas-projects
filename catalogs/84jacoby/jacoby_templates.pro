;+
; NAME:
;	JACOBY_TEMPLATES
;
; PURPOSE:
;       Generate templates to fit galaxy spectra based on the Jacoby
;       et al. 1984 stellar observations.
;
; CALLING SEQUENCE:
;
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;	
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMENTS:
;       We assume that the wavelengths are in the heliocentric frame
;       of reference and that the spectra have already been
;       de-reddened for interstellar extinction.  We remove the radial
;       velocity of each star based on measurements from the
;       literature.  
;
; EXAMPLE:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2002 January 28, U of A
;-

pro jacoby_templates, starflux, starwave, write=write, doplot=doplot

    light = 2.99792458D5 ; speed of light [km/s]

    eigenres = 4.5 ; Jacoby atlas resolution [Angstrom]
    
; read the jacoby atlas and sort by spectral type

    jhc = read_jacoby_atlas()

    srt = sort(strcompress(jhc.sp_type_lit,/remove))
    jhc = jhc[srt]
    
; only keep normal stars that have a radial velocity measurement 

    keep = where(strmatch(jhc.comment,'*normal star*') eq 1B or $
      strmatch(jhc.comment,'*cluster star*') eq 1B and $
      jhc.rv gt -900.0)
    jhc = jhc[keep]
    
; remove all O stars (lifetimes are too short), and GKM main-sequence
; stars (too faint)

    littype = strcompress(jhc.sp_type_lit,/remove)
    jhctype = strcompress(jhc.sp_type,/remove)

    keep = where((strmatch(littype,'O*') eq 0B) and (strmatch(littype,'G*V*') eq 0B) and $
      (strmatch(littype,'K*V*') eq 0B) and (strmatch(littype,'M*V*') eq 0B))
    jhc = jhc[keep]

; randomly remove stars of the same spectral type

    utype = uniq(littype[keep])
    jhc = jhc[utype]

; remove the following stars

    remstars = [$
      'HD12027',  $ ; pixel spike
      'HD23733',  $ ; stitching problem(s)
      'HD167451', $ ; pixel spike
      'HD14357',  $ ; small pixel spike
      'HD240344', $ ; strange near H-beta
      'HD13505',  $ ; stitching problem(s)
      'HD94028',  $ ; unwanted spectral type (F4V)
      'SAO102986',$ ; unwanted spectral type (F6V)
      'BD+051668',$ ; unwanted spectral type (M5V)
      'HD37767',  $ ; stitching problem(s)
      'HD39136',  $ ; stitching problem(s)
      'HD240296', $ ; stitching problem(s)
      'HD15866',  $ ; stitching problem(s)
      'SAO20899', $ ; stitching problem(s)
      'HD192832', $ ; stitching problem(s)
      'SAO87716', $ ; stitching problem(s)
      'HD842',    $ ; stitching problem(s)
      'HD186293', $ ; stitching problem(s)
      'HD1069',   $ ; stitching problem(s)
      'HD17647',  $ ; unwanted spectral type (G1V)
      'HD22193',  $ ; unwanted spectral type (G6V)
      'HD107214', $ ; unwanted spectral type (G2)
      'SAO57199', $ ; unwanted spectral type?? (F6V)
      'HD6111',   $ ; unwanted spectral type?? (F8V)
      'HD31084',  $ ; unwanted spectral type?? (F9V)
      'HD83140',  $ ; unwanted spectral type?? (F3IV)
      'HD260655', $ ; unwanted spectral type?? (M0V/K7)
      'BD+630137',$ ; unwanted spectral type?? (M1V/K5)
      'HD56030',  $ ; choose a F4III to represent this class (F6III)
      'SAO20603', $ ; choose a F4III to represent this class (F7III)
      'HD2506',   $ ; choose a G8III/G9III to represent this class (G4III)
      'HD112872', $ ; choose a G8III/G9III to represent this class (G6III)
      'SAO55164', $ ; choose a G8III/G9III to represent this class (K0III/G9III)
      'SAO37370', $ ; choose a A7I/F1II to represent this class (F0I/FOIb)
      'HD12842',  $ ; choose a A7I/F1II to represent this class (F3I/F3Ib)
      'HD17971'   $ ; choose a A7I/F1II to represent this class (F7I/F5Ia)
      ]

    doit = match_string(remstars,jhc.star,/exact,index=indx)
    keep = lindgen(n_elements(jhc))
    remove, indx, keep

    jhc = jhc[keep]

    srt = sort(jhc.jhc_id)
    jhc = jhc[srt]
    print_struct, jhc, ['jhc_id','star','sp_type','sp_type_lit','rv']

; ---------------------------------------------------------------------------

    starcount = n_elements(jhc)
    
    crval1 = float(jhc[0].wave[0])
    cd1_1 = float(jhc[0].wave[1]-jhc[0].wave[0])
    refwave = float(jhc[0].wave)
    npix = n_elements(refwave)

    starflux = fltarr(npix,starcount)
    
; initialize and fill the "output info" and the "plot info" structures 
          
    jhcinfo = {$
      template_name:         jhc.star,       $
      template_file:         jhc.star,       $
      template_type:         jhc.sp_type,    $ ; spectral type
      template_type_content: 'Spectral Type',$
      ntemplate:             starcount}
    
    for i = 0L, starcount-1L do begin

; correct for radial velocity
       
       flux = jhc[i].flux
       wave = jhc[i].wave
       rv = jhc[i].rv
       
       restwave = wave/sqrt((1+rv/light)/(1-rv/light)) ; rest wavelength
       starflux[*,i] = interpol(flux,restwave,refwave)
       
;      rest_cd1_1 = restwave[1]-restwave[0]
;      restflux = flux*cd1_1/rest_cd1_1
;      starflux[*,i] = interpol(restflux,restwave,refwave)
       
       if keyword_set(doplot) then begin
          djs_plot, refwave, starflux[*,i], xsty=3, ysty=3, charsize=2.0, charthick=2.0, $
            xtitle='Wavelength ('+angstrom()+')', ytitle='Relative f_{\lambda}', $
            xthick=2.0, ythick=2.0, yrange=[-0.1,1.1]
          legend, strcompress([string(jhc[i].jhc_id,format='(I0)'),jhc[i].star,$
            jhc[i].sp_type],/remove), /right, /top, box=0, $
            charsize=2.0, charthick=2.0
          cc = get_kbrd(1)
;          plot, wave, flux, ps=10, xr=[5000,5025], xsty=3, ysty=3
;          djs_oplot, restwave, newflux, ps=10, color='red'
;          djs_oplot, refwave, starflux[*,i], ps=10, color='green'
;          cc = get_kbrd(1)
       endif

    endfor 

    mkhdr, header, starflux
    sxaddpar, header, 'OBJECT', 'Jacoby stellar templates', im_today()
    sxaddpar, header, 'CRVAL1', crval1, ' central wavelength of first pixel'
    sxaddpar, header, 'CD1_1', cd1_1, ' dispersion [Angstrom/pixel]'
    sxaddpar, header, 'CRPIX1', 0, ' starting pixel (0-indexed)'
    sxaddpar, header, 'CTYPE1', 'LINEAR'
    sxaddpar, header, 'DC-FLAG', 0, ' log-linear flag'
    sxaddpar, header, 'EIGENRES', eigenres, 'resolution [Angstrom]'

; write out

    if keyword_set(write) then begin

       tpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
       tfile = 'jacoby_templates.fits'
       print, 'Writing '+tpath+tfile+'.'
       mwrfits, starflux, tpath+tfile, header, /create
       mwrfits, jhcinfo, tpath+tfile

    endif

return
end
