;+
; NAME:
;       KENN92_HTML_ATLAS
;
; PURPOSE:
;       Generate web page visualizations of the Kenn92 mergers. 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;       This routine assumes that all the appropriate directories
;       already exist.
;
; PROCEDURES USED:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004 May 04, U of A - written
;-

pro kenn92_html_atlas, object, write=write

; datapaths    
    
    webpath = atlas_path(/web)+'kenn92/'
    htmlpath = webpath+'html/'
    pngpath = webpath+'png/'

    baseref = '"http://cerebus.as.arizona.edu/~ioannis/kmatlas/kenn92/"'
    httphome = '"http://cerebus.as.arizona.edu/~ioannis"'
    httpstyle = '"http://cerebus.as.arizona.edu/~ioannis/styles/kmatlas.css"'
    
; restore all the fitting results

    splog, 'Reading the spectral fitting results.'
    atlas = read_kenn92(speclinefile=linefitfile,/silent)

    galaxy = strupcase(strcompress(atlas.galaxy,/remove))
    nicegalaxy = strupcase(atlas.nice_galaxy)
    nedgalaxy = strcompress(strupcase(atlas.ned_galaxy),/remove)
    ngalaxy = n_elements(galaxy)

; initialize individual galaxy file names
    
    htmlfile = 'html/'+strlowcase(galaxy)+'.html'
    asciiname = '../ascii/'+strlowcase(galaxy)+'.spec'
    fitsname = '../fits/'+strlowcase(galaxy)+'.fits'
    asciitablename = '../ascii/'+strlowcase(galaxy)+'.dat'
    fitstablename = '../fits/'+strlowcase(galaxy)+'_table.fits.gz'
    
    pngname = '../png/'+strlowcase(galaxy)+'.png'
    psname = '../png/'+strlowcase(galaxy)+'.ps.gz'

    specname = '../png/'+strlowcase(galaxy)+'_spec.png'
    specnameps = '../png/'+strlowcase(galaxy)+'_spec.ps.gz'

    spec1name = '../png/'+strlowcase(galaxy)+'_spec1.png'
    spec1nameps = '../png/'+strlowcase(galaxy)+'_spec1.ps.gz'
    spec2name = '../png/'+strlowcase(galaxy)+'_spec2.png'
    spec2nameps = '../png/'+strlowcase(galaxy)+'_spec2.ps.gz'
    spec3name = '../png/'+strlowcase(galaxy)+'_spec3.png'
    spec3nameps = '../png/'+strlowcase(galaxy)+'_spec3.ps.gz'
    spec4name = '../png/'+strlowcase(galaxy)+'_spec4.png'
    spec4nameps = '../png/'+strlowcase(galaxy)+'_spec4.ps.gz'

; ---------------------------------------------------------------------------
; KENN92ATLAS.HTML, KENN92TABLE.HTML, KENN92INTRO.HTML
; ---------------------------------------------------------------------------

    html_kenn92atlas, webpath, baseref, httphome, httpstyle, $
      nicegalaxy, htmlfile
    
; ---------------------------------------------------------------------------    
; individual galaxy web pages
; ---------------------------------------------------------------------------    
    
; basic properties
; ---------------------------------------------------------------------------

    bp = struct_trimtags(atlas,$
      select=['ra' ,'dec','cz'   ,'distance','rc3_type','pa','d25_maj','d25_min','inclination'],$
      format=['A15','A15','I10'  ,'F8.2'    ,'A10'     ,'I4','F7.2'   ,'F7.2'   ,'I4'         ])
    bplabels = ['RA','DEC','cz','d','T','PA','D25'    ,'d25'    ,'incl']
    bpunits = ['[J2000]','[J2000]','[km/s]','[Mpc]','&nbsp','[deg]','[arcmin]','[arcmin]','[deg]']
    nbptags = n_tags(bp)
    tagnames = tag_names(bp)
    
    for itag = 0L, nbptags-1L do begin
       prop = strtrim(reform(bp.(itag)),2)
       no = where(prop eq '',nno,comp=go,ncomp=ngo)
       if nno ne 0L then prop[no] = '-' else begin
          if (tagnames[itag] eq 'RC3_TYPE') then prop[go] = '-' else begin ; <-- NOT GENERAL!!!
             no = where(prop lt -900.0,nno,comp=go,ncomp=ngo)
             if nno ne 0L then prop[no] = '-'
             if ngo ne 0L then prop[go] = prop[go]
          endelse
       endelse
       bp.(itag) = strtrim(reform(prop,1,ngalaxy),2)
    endfor
    
; observing parameters
; ---------------------------------------------------------------------------

    op = struct_trimtags(atlas,$
      select=['driftdate','driftpa','driftap','driftscan','driftcode','driftphotflag','driftcomments'],$
      format=['A10','I3','I3','I3','I0','A1','A0'])
    oplabels = ['Date','Slit PA','Aperture','Scan','Code','Photometric?','Comments']
    opunits = ['&nbsp','[deg]','[arcsec]','[arcsec]','[123]','[Y/N]','&nbsp']
    noptags = n_tags(op)

; photometric properties
; ---------------------------------------------------------------------------

    pp = struct_trimtags(atlas,$
      select=['M_B' ,'U'   ,'B'   ,'V'   ,'J'   ,'H'   ,'Ks','IRAS_60','IRAS_100','FIR_LUM','fir_opt_ratio'],$
      format=['F7.1','F6.1','F6.1','F6.1','F6.1','F6.1','F6.1','F8.3','F8.3','F7.2','F8.2'])
    pplabels = ['M(B)','U','B','V','J','H','Ks','f(60)','f(100)','L(FIR)','L(FIR)/L(B)']
    ppunits = ['[mag]','[mag]','[mag]','[mag]','[mag]','[mag]','[mag]','[Jy]','[Jy]','[log L(Sun)]','&nbsp']
    npptags = n_tags(pp)

    for itag = 0L, npptags-1L do begin
       prop = strtrim(reform(pp.(itag)),2)
       no = where(prop lt -900.0,nno,comp=go,ncomp=ngo)
       if nno ne 0L then prop[no] = '-'
       if ngo ne 0L then prop[go] = prop[go]
       pp.(itag) = strtrim(reform(prop,1,ngalaxy),2)
    endfor

; emission-line properties
; ---------------------------------------------------------------------------

    linename = ['OII_3727','H_BETA','OIII_5007','OI_6300','H_ALPHA','NII_6584','SII_6716','SII_6731']
    linewavename = linename+'_WAVE'
    linechi2name = linename+'_CHI2'
    linenameew = linename+'_EW'
    linecontname = linename+'_CONTINUUM'
    linesigname = linename+'_SIGMA'

    linelabel = ['[O II]','H-BETA','[O III]','[O I]','H-ALPHA','[N II]','[S II]','[S II]']
    nline = n_elements(linelabel)

    spwave = struct_trimtags(atlas,select=linewavename,format=replicate('F7.2',nline))

    spchi2 = struct_trimtags(atlas,select=linechi2name,format=replicate('F12.1',nline))
    for itag = 0L, n_tags(spchi2)-1L do begin
       neg = where(spchi2.(itag) lt 0.0,nneg)
       if nneg ne 0L then spchi2[neg].(itag) = '-'
    endfor
    
; equivalent widths    
; ---------------------------------------------------------------------------
    
    spew = struct_trimtags(atlas,select=linenameew,format=replicate('F10.1',nline))
    for itag = 0L, n_tags(spew)-1L do begin
       flx = reform((spew.(itag))[0,*]) & err = reform((spew.(itag))[1,*]) & new = flx
       no = where(err eq -2.0,nno) & if (nno ne 0L) then new[no] = '-'
       no = where(err eq -3.0,nno) & if (nno ne 0L) then new[no] = '<'+strtrim(new[no],2)
       spew.(itag) = transpose([[new],[err]])
    endfor

; fluxes
; ---------------------------------------------------------------------------

    sptemp = struct_trimtags(atlas,select=linename)
    
    for itag = 0L, n_tags(sptemp)-1L do begin
       flx = reform((sptemp.(itag))[0,*]) & err = reform((sptemp.(itag))[1,*]) & new = flx
       go = where((err gt  0.0) or (err eq -3.0),ngo) & if (ngo ne 0L) then new[go] = alog10(new[go])
       sptemp.(itag) = transpose([[new],[err]])
    endfor

    spflx = struct_trimtags(sptemp,select=linename,format=replicate('F7.3',nline))
    for itag = 0L, n_tags(spflx)-1L do begin
       flx = reform((spflx.(itag))[0,*]) & err = reform((spflx.(itag))[1,*]) & new = flx
       no = where(err eq -2.0,nno) & if (nno ne 0L) then new[no] = '-'
       no = where(err eq -3.0,nno) & if (nno ne 0L) then new[no] = '<'+strtrim(new[no],2)
       spflx.(itag) = transpose([[new],[err]])
    endfor

; emission-line continuum
; ---------------------------------------------------------------------------

    sptemp = struct_trimtags(atlas,select=linecontname)
    
    for itag = 0L, n_tags(sptemp)-1L do begin
       flx = reform((sptemp.(itag))[0,*]) & err = reform((sptemp.(itag))[1,*]) & new = flx
       go = where((err gt  0.0) or (err eq -3.0),ngo) & if (ngo ne 0L) then new[go] = alog10(new[go])
       sptemp.(itag) = transpose([[new],[err]])
    endfor

    spcont = struct_trimtags(sptemp,select=linecontname,format=replicate('F7.3',nline))
    for itag = 0L, n_tags(spcont)-1L do begin
       flx = reform((spcont.(itag))[0,*]) & err = reform((spcont.(itag))[1,*]) & new = flx
       no = where(err eq -2.0,nno) & if (nno ne 0L) then new[no] = '-'
       no = where(err eq -3.0,nno) & if (nno ne 0L) then new[no] = '<'+strtrim(new[no],2)
       spcont.(itag) = transpose([[new],[err]])
    endfor

; emission-line sigma-width
; ---------------------------------------------------------------------------

    spsig = struct_trimtags(atlas,select=linesigname,format=replicate('F5.1',nline))

; stellar continuum properties
; ---------------------------------------------------------------------------

    schi2 = replicate({chi2: [0.0,0.0]},ngalaxy)
    schi2.chi2[0] = atlas.continuum_chi2
    schi2 = struct_trimtags(schi2,select='chi2',format='F10.3')
    
;   ebvc = replicate({starebv: [0.0,0.0]},ngalaxy)
;   ebvc.starebv[0,*] = reform(total(atlas.starebv,1)/atlas.ntemplate,1,ngalaxy)
;   ebvstr = struct_trimtags(ebvc,select='starebv',format='F10.3')
    
    cont = struct_trimtags(atlas,$
      select=['d4000_narrow','c41_50','lick_hb','lick_mgb','lick_fe','babs_h_alpha_ew',$
      'babs_h_beta_ew','babs_h_gamma_ew'],$
      format=['F10.3','F10.3','F10.3','F10.3','F10.3','F7.2','F7.2','F7.2'])
    cont = struct_addtags(schi2,cont)
;   cont = struct_addtags(ebvstr,cont)
    
    contlabels = ['Chi2','Dn(4000)','41-50','Lick Hb','Lick Mg b','Lick <Fe>','EW(Ha) Abs','EW(Hb) Abs','EW(Hg) Abs']
    contunits = ['&nbsp','&nbsp','[mag]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]']
;   contlabels = ['E(B-V)','Dn(4000)','41-50','Lick Hb','Lick Mgb','Lick Fe','EW(Ha) Abs','EW(Hb) Abs','EW(Hg) Abs']
;   contunits = ['[mag]','&nbsp','[mag]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]','[Angstrom]']
    nconttags = n_tags(cont)

    for itag = 0L, nconttags-1L do begin
       flx = reform((cont.(itag))[0,*]) & err = reform((cont.(itag))[1,*]) & new = flx
       no = where(err eq -1.0,nno) & if (nno ne 0L) then new[no] = '-'
       cont.(itag) = transpose([[new],[err]])
    endfor

    ntemp = atlas[0].ntemplate
    temp = struct_trimtags(atlas,select=['template_age','continuum_fraction_V',$
      'continuum_mass_fraction','continuum_ebv'])
    temp.(0) = temp.(0)/1E6
    temp.(1) = temp.(1)*100.0
    temp.(2) = temp.(2)*100.0
    temp = struct_trimtags(temp,select=tag_names(temp),format=['I7','I3','I3','F10.3'])

    for i = 0L, ngalaxy-1L do begin
       junk = temp[i].(1)
       zero = where(junk eq 0.0,nzero)
       if nzero ne 0L then begin
          junk[zero] = '-'
          temp[i].(1) = junk
       endif
       junk = temp[i].(2)
       zero = where(junk eq 0.0,nzero)
       if nzero ne 0L then begin
          junk[zero] = '-'
          temp[i].(2) = junk
       endif
    endfor
    
; ---------------------------------------------------------------------------
; write out
; ---------------------------------------------------------------------------
    
    pushd, webpath
    
    splog, 'Writing individual galaxy web pages.'
    for i = 0, ngalaxy-1L do begin

       openw, lun, htmlfile[i], /get_lun
       printf, lun, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional'+$
         '//EN" '+httphome+'>' 
       printf, lun, '<html>'
       printf, lun, '<head>'
       printf, lun, '<link rel="stylesheet" type="text/css" href='+httpstyle+' />'
       printf, lun, '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></meta>'
;      printf, lun, '<title>'+nicegalaxy[i]+'</title>'
       printf, lun, '</head>'
       printf, lun, '<body>'
       printf, lun, '<h1><a href="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?objname='+nedgalaxy[i]+$
         '&extend=no&out_csys=Equatorial&out_equinox=B2000.0&obj_sort=RA+or+Longitude&zv=z&zv_breaker=10000.0" target="_top">'+$
         nicegalaxy[i]+'</a></h1>'

       printf, lun, '<div class="menu right">'
       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
       printf, lun, '<a href="'+backfile+'">Back</a>'
       printf, lun, '</div>'

       printf, lun, '<div class="menu left">'
       printf, lun, '<a href="../../index.html" target="_top">Home</a><br />'
       printf, lun, '<a href="../kenn92intro.html" target="galaxy">Intro</a>'
       printf, lun, '</div><br /><br clear="bottom"/>'

; basic properties table       
; ---------------------------------------------------------------------------
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Basic Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, nbptags-1L do printf, lun, '<th>'+bplabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, nbptags-1L do printf, lun, '<th><span class="smaller">'+bpunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, nbptags-1L do printf, lun, '<td>'+bp[i].(itag)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br />'

; image and spectrum
; ---------------------------------------------------------------------------

       printf, lun, '<table width="100%" border="0" cellspacing="1" cellpadding="1">'
       printf, lun, '<tbody>'
       printf, lun, '<tr><td width="49%" rowspan="2" valign="top"><a style="color: silver;" href="'+$
         psname[i]+'"><img width="99%" src="'+pngname[i]+'" alt="'+galaxy[i]+' Thumbnail" /></a></td>'
       printf, lun, '<td width="49%" valign="top"><a style="color: silver;" href="'+specnameps[i]+'"><img width="99%" src="'+$
         specname[i]+'" alt="'+galaxy[i]+' Thumbnail" /></a></td></tr>'
       printf, lun, '<tr><td width="100%" align="center" valign="center">'
       printf, lun, '   <table style="font-size: 150%;" border="1">'
       printf, lun, '   <tr>'
       printf, lun, '      <td>Observed Spectrum</td><td><a href="'+spec1name[i]+'">png</a> <a href="'+$
         spec1nameps[i]+'">ps</a></td>'
       printf, lun, '      <td>Continuum Model</td><td><a href="'+spec2name[i]+'">png</a> <a href="'+$
         spec2nameps[i]+'">ps</a></td>'
       printf, lun, '   </tr>'
       printf, lun, '   <tr>'
       printf, lun, '      <td>Emission Line Model</td><td><a href="'+spec3name[i]+'">png</a> <a href="'+$
         spec3nameps[i]+'">ps</a></td>'
       printf, lun, '      <td>Composite Model</td><td><a href="'+spec4name[i]+'">png</a> <a href="'+$
         spec4nameps[i]+'">ps</a></td>'
       printf, lun, '   </tr>'
       printf, lun, '   </table>'
       printf, lun, '</td></tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br clear="bottom"/>'

; observations table
; ---------------------------------------------------------------------------

;      printf, lun, '<div class="tablecenter">'
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Observing Parameters</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, noptags-1L do printf, lun, '<th>'+oplabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, noptags-1L do printf, lun, '<th><span class="smaller">'+opunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, noptags-1L do printf, lun, '<td>'+op[i].(itag)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
;      printf, lun, '</div>'
       printf, lun, '<br /><br />'

; photometric properties table       
; ---------------------------------------------------------------------------
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Photometric Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<th>'+pplabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<th><span class="smaller">'+ppunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, npptags-1L do printf, lun, '<td>'+pp[i].(itag)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br /><br />'

; emission-line properties table       
; ---------------------------------------------------------------------------
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Emission Line Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       printf, lun, '<th>&nbsp</th>' ; blank
       for iline = 0L, nline-1L do printf, lun, '<th>'+linelabel[iline]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">Wavelength<br /><span class="smaller">[Angstrom]</span></td>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spwave[i].(iline)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">EW<br /><span class="smaller">[Angstrom]</span></th>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spew[i].(iline)[0,*]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">Flux<br /><span class="smaller">[log erg/s/cm2]</span></th>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spflx[i].(iline)[0,*]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">Continuum<br /><span class="smaller">[log erg/s/cm2/A]</span></th>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spcont[i].(iline)[0,*]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">Gaussian Width<br /><span class="smaller">[km/s]</span></th>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spsig[i].(iline)[0,*]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       printf, lun, '<th class="rowhead">Reduced Chi2</td>'
       for iline = 0L, nline-1L do printf, lun, '<td>'+spchi2[i].(iline)+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br /><br />'

; continuum properties table       
; ---------------------------------------------------------------------------
       
       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<caption>Continuum Properties</caption>'
       printf, lun, '<thead>'
       printf, lun, '<tr>'
       for itag = 0L, nconttags-1L do printf, lun, '<th>'+contlabels[itag]+'</th>'
       printf, lun, '</tr>'
       printf, lun, '<tr>'
       for itag = 0L, nconttags-1L do printf, lun, '<th><span class="smaller">'+contunits[itag]+'</span></th>'
       printf, lun, '</tr>'
       printf, lun, '</thead>'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       for itag = 0L, nconttags-1L do printf, lun, '<td>'+(cont[i].(itag))[0]+'</td>'
       printf, lun, '</tr>'
       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br />'

       printf, lun, '<table width="100%" border="2" cellspacing="0" cellpadding="2" class="props">'
       printf, lun, '<tbody>'
       printf, lun, '<tr>'
       printf, lun, '<th>Age <span class="smaller">[Myr]</span></th>'
       for itemp = 0L, ntemp-1L do printf, lun, '<td>'+(temp[i].(0))[itemp]+'</td>'
       printf, lun, '</tr>'

       printf, lun, '<tr>'
       printf, lun, '<th>Light Fraction (V-band) <span class="smaller">[%]</span></th>'
       for itemp = 0L, ntemp-1L do printf, lun, '<td>'+(temp[i].(1))[itemp]+'</td>'
       printf, lun, '</tr>'

       printf, lun, '<tr>'
       printf, lun, '<th>Mass Fraction <span class="smaller">[%]</span></th>'
       for itemp = 0L, ntemp-1L do printf, lun, '<td>'+(temp[i].(2))[itemp]+'</td>'
       printf, lun, '</tr>'

       printf, lun, '<tr>'
       printf, lun, '<th>E(B-V) <span class="smaller">[mag]</span></th>'
       for itemp = 0L, ntemp-1L do printf, lun, '<td>'+(temp[i].(3))[itemp]+'</td>'
       printf, lun, '</tr>'

       printf, lun, '</tbody>'
       printf, lun, '</table>'
       printf, lun, '<br /><br />'

;;; download       
;;; ---------------------------------------------------------------------------
;;       
;;       printf, lun, '<div class="box download">'
;;       printf, lun, '<h2>Download</h2>'
;;       printf, lun, '<ul>'
;;       printf, lun, '<li><a href="'+fitsname[i]+'">FITS  Spectrum</a></li>'
;;       printf, lun, '<li><a href="'+asciiname[i]+'">ASCII Spectrum</a></li>'
;;       printf, lun, '<li><a href="'+fitstablename[i]+'">FITS Data Table</a></li>'
;;       printf, lun, '</ul>'
;;       printf, lun, '</div>'
;;       
;;       printf, lun, '<p class="menu right">'
;;       if i eq ngalaxy-1L then nextfile = '../'+htmlfile[0] else nextfile = '../'+htmlfile[i+1L]
;;       printf, lun, '<a href="'+nextfile+'">Next</a><br />'
;;       if i eq 0L then backfile = '../'+htmlfile[ngalaxy-1L] else backfile = '../'+htmlfile[i-1L]
;;       printf, lun, '<a href="'+backfile+'">Back</a>'
;;       printf, lun, '</p>'
;; 
;;       printf, lun, '<p class="menu left">'
;;       printf, lun, '<a href="../../index.html" target="_top">Home</a><br />'
;;       printf, lun, '<a href="../kenn92intro.html" target="galaxy">Intro</a>'
;;       printf, lun, '</p><br /><br clear="bottom"/><br />'

       printf, lun, '<p class="email">Last update '+im_today()+'.  Email <a href="mailto:jmoustakas@as.arizona.edu">'+$
         'J. Moustakas</a> with questions or comments.</p>'

       printf, lun, '</body>'
       printf, lun, '</html>'
       free_lun, lun
       
    endfor 

    popd
    
return
end    
