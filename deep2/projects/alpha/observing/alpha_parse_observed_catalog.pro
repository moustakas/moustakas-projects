pro alpha_parse_observed_catalog
; jm08apr16nyu - 
; jm08jul12nyu - read the catalog of observed objects, cross-match
; with the various parent catalogs, and write out some useful info on
; the properties of each object; this has to be done this way because
; I didn't quite version-control the parent catalogs of objects that
; were observed in nov07 vs apr08; plus some objects were added
; on-the-fly at the telescope
; jm09jan09nyu - updated to use the sep08.cat
    
    common sdss, sdsscat, sdssspec, sdss_mpa
    
    latest = 'sep08'
    path = deep2_path(/projects)+'alpha/'
    
; read the latest catalog

    readcol, path+latest+'.cat', parent_id, galaxy, ra, $
      dec, format='I,A,A,A', /silent
    nobj = n_elements(parent_id)

    cat = {$
      parent_id: -1L, $
      galaxy:    '', $
      ra:      0.0D, $
      dec:     0.0D, $
      z:       -1.0, $
      mjd:      -1L, $
      plate:     -1, $
      fiber:     -1}
;     mg:      -1.0, $
;     gr:      -1.0}
    cat = replicate(cat,nobj)

    cat.parent_id = parent_id
    cat.galaxy = strtrim(galaxy,2)
    cat.ra = im_hms2dec(ra)*15.0D
    cat.dec = im_hms2dec(dec)

; cross-match the catalog with itself to find duplicates: obj#1 and
; obj#79 are duplicates
;   ss = where(strmatch(galaxy,'*sdss*',/fold))
;   ingroup = spheregroup(cat[ss].ra,cat[ss].dec,2.0/3600.0,$
;     multgroup=mult,first=first,next=next)
    
; get the SDSS redshift
    if (n_elements(sdsscat) eq 0L) then sdsscat = hogg_mrdfits($
      vagc_name('object_sdss_imaging'),1,nrow=28800L);,range=[0,10])
    if (n_elements(sdssspec) eq 0L) then sdssspec = hogg_mrdfits($
      vagc_name('object_sdss_spectro'),1,nrow=28800L);,range=[0,10])

; spherematch with maxmatch=0 because there is one duplicate    
    spherematch, sdsscat.ra, sdsscat.dec, cat.ra, cat.dec, $
      3.0D/3600.0, m1, m2, maxmatch=0
    cat[m2].z = sdssspec[m1].z
    cat[m2].mjd = sdssspec[m1].mjd
    cat[m2].plate = sdssspec[m1].plate
    cat[m2].fiber = sdssspec[m1].fiberid

;   need = where(strmatch(galaxy,'*sdss*',/fold) and cat.z lt 0.0,nneed)
;   struct_print, cat[need]

; total hack!    
    deep1 = mrdfits(path+'deep2_alpha_08apr.fits.gz',1)
    deep2 = mrdfits(path+'deep2_alpha_08sep.fits.gz',1)
    deep = [struct_trimtags(deep1,sel=['galaxy','z','ra','dec']),$
      struct_trimtags(deep2,sel=['galaxy','z','ra','dec'])]
;   if (n_elements(deep) eq 0L) then deep = mrdfits(deep2_path(/analysis)+$
;     'deep2_kcorr_'+deep2_version(/kcorr)+'.fits.gz',1)

; match against DEEP2
    spherematch, cat.ra, cat.dec, deep.ra, deep.dec, $
      5.0D/3600.0, m1_deep, m2_deep
    cat[m1_deep].z = deep[m2_deep].z
;   cat[m1_deep].mg = deep[m2_deep].ugriz_absmag[1]
;   cat[m1_deep].gr = deep[m2_deep].ugriz_absmag[1]-deep[m2_deep].ugriz_absmag[2]
    splog, 'Found '+string(n_elements(m1_deep),format='(I0)')+$
      ' targets from DEEP.'

    need = where((strmatch(galaxy,'*star*',/fold) eq 0B) and cat.z lt 0.0,nneed)
    if (nneed ne 0L) then begin
       splog, 'Oh no!'
       struct_print, cat[need]
    endif

    struct_print, cat
    mwrfits, cat, path+'alpha_observed_catalog_'+latest+'.fits', /create

stop    
    
return
end

;    sdss1 = mrdfits('~/home/research/data/sdss/mpa_dr7_v5_2/gal_info_dr7_v5_2.fit.gz',1)
;    ww = where(sdss1.ra gt -9000.0 and sdss1.dec gt -9000.0)
;    sdss = sdss1[ww]
;
;    spherematch, sdss.ra, sdss.dec, cat.ra, cat.dec, 5.0D/3600.0, m1, m2, maxmatch=1
;    cat[m2].z = sdss[m1].z
;
;    need = where(strmatch(galaxy,'*sdss*',/fold) and cat.z lt 0.0,nneed)
;    struct_print, cat[need]
;
;    djs_plot, sdss.ra, sdss.dec, ps=3
;    djs_oplot, cat[ss].ra, cat[ss].dec, ps=4, color='red'
;    plots, cat[need].ra, cat[need].dec, ps=7, sym=2, color=djs_icolor('dark green')
;    
;    if (n_elements(sdss_mpa) eq 0L) then sdss_mpa = read_sdss_mpa(/info)
;
;; match against the full SDSS/MPA catalog
;    need = where(strmatch(galaxy,'*sdss*',/fold) and cat.z lt 0.0,nneed)
;    struct_print, cat[need]
;    if (nneed ne 0L) then begin
;       spherematch, sdss_mpa.ra, sdss_mpa.dec, cat[need].ra, cat[need].dec, $
;         5.0D/3600.0, m1, m2
;       cat[need[m2]].z = sdss_mpa[m1].z
;       splog, 'Found '+string(n_elements(m1_mpa),format='(I0)')+$
;         ' targets from SDSS/MPA.'
;    endif
;    
;    need = where((strmatch(galaxy,'*star*',/fold) eq 0B) and cat.z lt 0.0,nneed)
;    struct_print, cat[need]
;
;; get the redshift - fix this!!    
;    for cc = 0L, n_elements(cat)-1L do begin
;       name = strsplit(galaxy[cc],'.',/extract)
;       if (strtrim(name[0],2) eq 'sdss') then begin
;          cat[cc].mjd = name[1]
;          cat[cc].plate = name[2]
;          cat[cc].fiber = name[3]
;       endif
;    endfor
;
;    ww = where(cat.plate gt 0)
;    readspec, cat[ww].plate, cat[ww].fiber, mjd=cat[ww].mjd, zans=zz
;    cat[ww].z = zz.z
;    niceprint, cat[ww].ra, zz.plug_ra, cat[ww].dec, zz.plug_dec, cat[ww].mjd, zz.mjd
;
;; this is a complete hack!    
;;   sdss1 = mrdfits(path+'sdss_alpha_07nov.fits.gz',1)
;;   sdss2 = mrdfits(path+'sdss_alpha_08apr.fits.gz',1)
;;   sdss3 = mrdfits(path+'sdss_alpha_08apr_dr6.fits.gz',1)
;;   sdss4 = mrdfits(path+'sdss_alpha_08apr_crap.fits.gz',1)
;;   sdss_dr6 = [struct_trimtags(sdss1,sel=['galaxy','z','ra','dec']),$
;;     struct_trimtags(sdss2,sel=['galaxy','z','ra','dec']),$
;;     struct_trimtags(sdss3,sel=['galaxy','z','ra','dec']),$
;;     struct_trimtags(sdss4,sel=['galaxy','z','ra','dec'])]
;;   if (n_elements(sdss_dr6) eq 0L) then sdss_dr6 = $
;;     hogg_mrdfits('/mount/moon1/ioannis/sdss_dr4_v5_1b/gal_info_dr6_v5_1.fit',1)
;;   if (n_elements(sdss_dr6) eq 0L) then sdss_dr6 = $
;;     hogg_mrdfits('/mount/moon1/ioannis/sdss_dr4_v5_1b/sdss_main_kcorr_dr6.fits',1)
;;   if (n_elements(sdss_dr6) eq 0L) then sdss_dr6 = read_sdss_main(/kcorr,/dr6)
;;   if (n_elements(sdss_dr6) eq 0L) then sdss_dr6 = read_sdss_vagc_mpa($
;;     sample='dr6',poststr='0',/kcorr)
;
;; match against the full SDSS/DR6 catalog
;;   spherematch, cat.ra, cat.dec, sdss_dr6.ra, sdss_dr6.dec, $
;;     5.0D/3600.0, m1_dr6, m2_dr6
;;   cat[m1_dr6].z = sdss_dr6[m2_dr6].z
;;   cat[m1_dr6].mg = sdss_dr6[m2_dr6].ugriz_absmag[1]
;;   cat[m1_dr6].gr = sdss_dr6[m2_dr6].ugriz_absmag[1]-sdss_dr6[m2_dr6].ugriz_absmag[2]
;;   splog, 'Found '+string(n_elements(m1_dr6),format='(I0)')+$
;;     ' targets from SDSS/DR6.'
;    
;;   need = where((strmatch(galaxy,'*star*',/fold) eq 0B) and cat.z lt 0.0,nneed)
;;   struct_print, cat[need]
    
