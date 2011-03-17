pro ages_arjun_comparison, ages, arjun, agesdust, ancillary, write=write
; jm06feb23uofa - 
    
    datapath = ages_path(/mass)

    if (n_elements(ages) eq 0L) then ages = ages_read_ancillary()
    if (n_elements(arjun) eq 0L) then arjun = rsex(datapath+'Redgalaxies_with_redshifts.dat')

;   zmerge = mrdfits(ages_path(/analysis)+'catalog.zmerge.fits.gz',1,/silent)
;   zindx = cmset_op(zmerge.zmerge_catid,'and',arjun.ages_index,/index)

    match, ages.ages_id, arjun.ages_index, indx1, indx2
    fitages = ages[indx1] & fitarjun = arjun[indx2]
;   help, fitarjun, arjun
;   plot, fitages.z, fitarjun.zspec, ps=4

    zobj = fitages.z
    galaxy = strtrim(fitages.galaxy,2)
    ngalaxy = n_elements(fitages)

;   missing = lindgen(n_elements(arjun)) & remove, indx2, missing
;   missing_index = arjun[missing].ages_index

; ---------------------------------------------------------------------------
; write out D4000
; ---------------------------------------------------------------------------

    if (n_elements(agesdust) eq 0L) then agesdust = read_ages()
    fitagesdust = agesdust[indx1]
    if (n_elements(ancillary) eq 0L) then ancillary = ages_read_ancillary()
    fitancillary = ancillary[indx1]
;   outarjun = struct_addtags(fitarjun,struct_trimtags(agesdust[indx1],$
;     select=['Z_OBJ','Z_LINE','D4000','D4000_MODEL','D4000_NARROW','D4000_NARROW_MODEL']))
;   mwrfits, outarjun, datapath+'ages_redgalaxies_d4000.fits', /create
;;  openw, lun, datapath+'ages_redgalaxies_d4000.dat', /get_lun
;;  struct_print, outarjun, lun=lun
;;  free_lun, lun

; ---------------------------------------------------------------------------
; photometric bands of interest
; ---------------------------------------------------------------------------

    modelsprefix = 'ages_BwRIKJKs'
    massprefix = 'ages_Redgalaxies_BwRIKJKs'
    filterlist = ['dwfs_Bw','dwfs_R','dwfs_I','dwfs_onis_K','flamex_flamingos_J','flamex_flamingos_Ks']+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilter = n_elements(filterlist)
    
; ---------------------------------------------------------------------------
; convert the photometry in each band to MAGGIES and do the fit!
; ---------------------------------------------------------------------------

    bands = 'phot_'+['bw','r','i','k','flamj','flamk']
    maggies = convert_to_maggies(ages,bands,filtinfo.vega2ab,maggies_invvar=maggies_invvar)

    fitindx = lindgen(ngalaxy)
    thesefilters = filterlist

;   isedfit,maggies[*,fitindx],maggies_invvar[*,fitindx],zobj[fitindx],result,filterlist=thesefilters,$
;     modelsprefix=modelsprefix,massprefix=massprefix,galaxy=galaxy[fitindx],/write
    isedfit_qaplot, result, result_info, datapath=datapath, massprefix=massprefix, /postscript

;   isedfit_measure, measure, massprefix=massprefix, /write

    w = where(measure.d4000_narrow gt 0.0)
    plot, measure[w].d4000_narrow, fitagesdust[w].d4000_narrow_model[0], ps=4, xr=[0.8,2.5], yr=[0.8,2.5]

;   ages_display_spectrum, fitagesdust[0:10], fitancillary[0:10], specfit=specfit
    
stop    
    
return
end
    
; ---------------------------------------------------------------------------    

;  for j=0,n_elements(result)-1 do for i=0,n_tags(result)-1 do if size(result[j].(i),/type) ne 7L then if total(finite(result[j].(i)) eq 0B) gt 0.0 then print, j, i

;;    datapath = ages_path(/before_skysubpca)
;;    outpath = ages_path(/after_skysubpca)
;;    analysis_path = ages_path(/analysis)
;;
;;    cat = mrdfits(analysis_path+'catalog.cat.noguidestars.fits.gz',1,/silent)
;;
;;    cat_arjun = cat[arjun.ages_index]
;;    
;;    datafiles = file_search(datapath+'spectra_???.fits.gz',count=npass)
;;    infofiles = file_search(datapath+'target_info_???.fits.gz',count=ninfo)
;;    for ipass=0L, npass-1L do if (n_elements(allinfo) eq 0L) then $
;;      allinfo = mrdfits(infofiles[ipass],1,/silent) else $
;;      allinfo = struct_append(allinfo,mrdfits(infofiles[ipass],1,/silent))
;;    help, allinfo
;;
;;    searchrad = 1.5
;;    ntot = im_djs_angle_match(cat.cat_ra*15.0,cat.cat_dec,$
;;;   ntot = im_djs_angle_match(allinfo.zmerge_ra,allinfo.zmerge_dec,$
;;      arjun.ra,arjun.dec,dtheta=searchrad/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist,mmax=mmax)
;;    match = where(mindx ne -1L,comp=nomatch)
;;    splog, 'Found '+string(ntot,format='(I0)')+' matching objects.'
;;
;;    missingindx = lindgen(n_elements(arjun))
;;    match, allinfo.zmerge_catid, arjun.ages_index, indx1, indx2
;;;   niceprint, allinfo[indx1].z, arjun[indx2].zspec
;;    remove, indx2, missingindx
;;    help, arjun, indx1, missingindx
;;    struct_print, arjun[missingindx]

; ---------------------------------------------------------------------------    
