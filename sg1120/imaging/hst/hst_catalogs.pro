pro hst_catalogs
; RELEGATED!!

       assoc_radius = 0.5      ; [arcsec]
       assoc_radius_zcat = 1.0 ; [arcsec]

       catname_hst = catpath+'sg1120_hst.cat'
;      catname_hst = catpath+'sg1120_hst_merged.cat'
       catname_gprime = catpath+'sg1120_gprime.cat'
       catname_rprime = catpath+'sg1120_rprime.cat'
       catname_zcat = catpath+'sg1120_07jun25.zcat.sex'

       outcatname_hst    = repstr(catname_hst,'_hst_','_hst_matched_')
       outcatname_gprime = repstr(catname_gprime,'_gprime_','_gprime_matched_')
       outcatname_rprime = repstr(catname_rprime,'_rprime_','_rprime_matched_')

       outcatname_hst_zcat    = repstr(catname_hst,'_hst_','_hst_zmatched_')
       outcatname_gprime_zcat = repstr(catname_gprime,'_gprime_','_gprime_zmatched_')
       outcatname_rprime_zcat = repstr(catname_rprime,'_rprime_','_rprime_zmatched_')

       if (file_test(catname_hst,/regular) eq 0L) then begin
          splog, 'HST catalog '+catname_hst+' not found.'
          return
       endif
       splog, 'Reading '+catname_hst
       hstcat = rsex(catname_hst)

       if (file_test(catname_rprime,/regular) eq 0L) then begin
          splog, 'r-band LDSS3 catalog '+catname_rprime+' not found.'
          return
       endif
       rprime = rsex(catname_rprime)

       if (file_test(catname_gprime,/regular) eq 0L) then begin
          splog, 'g-band LDSS3 catalog '+catname_gprime+' not found.'
          return
       endif
       gprime = rsex(catname_gprime)

       if (file_test(catname_zcat,/regular) eq 0L) then begin
          splog, 'Redshift catalog '+catname_zcat+' not found.'
          return
       endif
       splog, 'Reading '+catname_zcat
       zcat = rsex(catname_zcat)

; match the HST catalog to the r-prime catalog; the output should be
; identical to matching to the g-prime catalog, since they were
; generated in double-image mode, detecting on the r-band image       
       
       splog, 'Matching HST catalog to the r-prime catalog.'
       spherematch, hstcat.xwin_world, hstcat.ywin_world, rprime.xwin_world, $
         rprime.ywin_world, assoc_radius/3600.0, hstmatch, grmatch, dist12, $
         maxmatch=1

       splog, 'Writing '+outcatname_hst
       wsex, hstcat[hstmatch], outfile=outcatname_hst
       splog, 'Writing '+outcatname_gprime
       wsex, gprime[grmatch], outfile=outcatname_gprime
       splog, 'Writing '+outcatname_rprime
       wsex, rprime[grmatch], outfile=outcatname_rprime

; also write out a stellar catalog

;      plot, hstcat[hstmatch].mag_auto, hstcat[hstmatch].class_star, ps=4, xsty=3, ysty=3, yr=[0,1]
;      oplot, 22.0*[1,1], !y.crange & oplot, !x.crange, 0.8*[1,1]
;      plot, rprime[grmatch].mag_auto, rprime[grmatch].class_star, ps=4, xsty=3, ysty=3, xr=[-15,0], yr=[0,1]
;      oplot, -9.5*[1,1], !y.crange & oplot, !x.crange, 0.8*[1,1]
;      plot, gprime[grmatch].mag_auto, gprime[grmatch].class_star, ps=4, xsty=3, ysty=3, xr=[-15,0], yr=[0,1]
;      oplot, -9.5*[1,1], !y.crange & oplot, !x.crange, 0.8*[1,1]

       stars = where((hstcat[hstmatch].class_star gt 0.8) and (hstcat[hstmatch].mag_auto lt 21.0),nstars)
;      stars = where((rprime[grmatch].class_star gt 0.8) and (rprime[grmatch].mag_auto gt -9.5),nstars)

;      plothist, rprime[grmatch[stars]].fwhm_image*0.188, bin=0.05, xr=[0,1.5]
;      plothist, gprime[grmatch[stars]].fwhm_image*0.188, bin=0.05, xr=[0,1.5]
;      stats = im_stats(rprime[grmatch[stars]].fwhm_image*0.188,/verbose,sigrej=3.0)
;      stats = im_stats(gprime[grmatch[stars]].fwhm_image*0.188,/verbose,sigrej=3.0)
;      struct_print, struct_trimtags(rprime[grmatch[stars]],select=['X_IMAGE','Y_IMAGE']), /no_head, file='junk.reg'
;      struct_print, struct_trimtags(hstcat[hstmatch[stars]],select=['X_IMAGE','Y_IMAGE']), /no_head, file='junk2.reg'

; --------------------------------------------------       
;      splog, 'Writing '+repstr(outcatname_hst,'.cat','_STARS.cat')
;      wsex, hstcat[hstmatch[stars]], outfile=repstr(outcatname_hst,'.cat','_STARS.cat')
;      splog, 'Writing '+repstr(outcatname_gprime,'.cat','_STARS.cat')
;      wsex, gprime[grmatch[stars]], outfile=repstr(outcatname_gprime,'.cat','_STARS.cat')
;      splog, 'Writing '+repstr(outcatname_rprime,'.cat','_STARS.cat')
;      wsex, rprime[grmatch[stars]], outfile=repstr(outcatname_rprime,'.cat','_STARS.cat')
; --------------------------------------------------       

; match the HST-matched catalogs to the redshift catalog
       
       splog, 'Matching the HST-matched catalogs to the redshift catalog.'
       spherematch, hstcat[hstmatch].xwin_world, hstcat[hstmatch].ywin_world, $
         zcat.ra, zcat.dec, assoc_radius_zcat/3600.0, hstmatch_zcat, zcatmatch, $
         dist12, maxmatch=1

       zcat_out = im_struct_trimtags(zcat[zcatmatch],select=tag_names(zcat[0]),$
         newtags='ZCAT_'+tag_names(zcat[0]))

       splog, 'Writing '+outcatname_hst_zcat
       wsex, struct_addtags(hstcat[hstmatch[hstmatch_zcat]],zcat_out), outfile=outcatname_hst_zcat
       splog, 'Writing '+outcatname_gprime_zcat
       wsex, struct_addtags(gprime[grmatch[hstmatch_zcat]],zcat_out), outfile=outcatname_gprime_zcat
       splog, 'Writing '+outcatname_rprime_zcat
       wsex, struct_addtags(rprime[grmatch[hstmatch_zcat]],zcat_out), outfile=outcatname_rprime_zcat
    
    
; ----------------------------------------------------------------------------
; generate the master HST catalog
; ----------------------------------------------------------------------------
    
    if keyword_set(hst) then begin

       imagelist = mosaicpath+'sg1120_hst.fits' ; 0.035 arcsec/pixel
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catname = catpath+'sg1120_hst_'+date+'.cat'

       pixscale = 0.03 ; 0.035
       photaper = [1.0,2.0,3.0] ; [diameter,arcsec]
       photaper = strjoin(strtrim(string(photaper/pixscale,format='(F10.2)'),2),',') ; [diameter,pixel]
       seeing = '0.105'
       magzero = string(25.501+2.5*alog10(2000.0),format='(F7.4)')

       checkimage_type = 'NONE' ; 'SEGMENTATION'
       checkimage_name = mosaicpath+'hst.segmentation.fits'

       spawn, 'sex '+imagelist+' -c '+sexconfig+' -PARAMETERS_NAME '+sexparam+$
         ' -FILTER_NAME '+sexconv_hst+' -STARNNW_NAME '+sexnnw+' -CATALOG_NAME '+catname+' -CATALOG_TYPE ASCII_HEAD'+$
         ' -DETECT_THRESH 1.5 -DETECT_MINAREA 9 -DEBLEND_NTHRESH 16 -DEBLEND_MINCONT 0.015'+$
         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist+' -WEIGHT_THRESH 0 -WEIGHT_GAIN N'+$
         ' -PHOT_APERTURES '+photaper+' -MAG_ZEROPOINT '+magzero+' -SEEING_FWHM '+seeing+$
         ' -BACK_SIZE 512 -BACKPHOTO_TYPE LOCAL -BACKPHOTO_THICK 26 -INTERP_TYPE NONE -GAIN 2.0'+$
         ' -VERBOSE_TYPE NORMAL -NTHREADS 4 '+$
;        ' -MEMORY_OBJSTACK 65536 -MEMORY_PIXSTACK 10000000 -MEMORY_BUFSIZE 65534 -VERBOSE_TYPE NORMAL -NTHREADS 2 '+$
         ' -CHECKIMAGE_TYPE '+checkimage_type+' -CHECKIMAGE_NAME '+checkimage_name, /sh

    endif

stop    
    
return
end
    
;      spawn, 'sex '+strjoin(imagelist,',')+' -c '+sexconfig+' -PARAMETERS_NAME '+sexparam+$
;        ' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
;        ' -CATALOG_NAME '+catname+' -CATALOG_TYPE ASCII_HEAD -DETECT_MINAREA 5'+$
;        ' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER Y -FILTER_THRESH " "'+$
;        ' -DEBLEND_NTHRESH 32 -DEBLEND_MINCONT 0.005 -CLEAN Y -CLEAN_PARAM 1.0 -MASK_TYPE CORRECT'+$
;        ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+strjoin(weightlist,',')+' -WEIGHT_THRESH 0'+$
;        ' -PHOT_APERTURES '+photaper+' -MAG_ZEROPOINT '+magzero+' -SEEING_FWHM '+seeing+$
;        ' -BACK_TYPE AUTO -BACK_SIZE 64 -BACK_FILTERSIZE 3 -BACKPHOTO_TYPE GLOBAL '+$
;        ' -BACKPHOTO_THICK 24 -BACK_FILTTHRESH 0.0 -INTERP_TYPE NONE '+$
;        ' -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 2 '+$
;        ' -CHECKIMAGE_TYPE '+checkimage_type+' -CHECKIMAGE_NAME '+checkimage_name, /sh

;; generate the LDSS3 catalogs using the master/merged HST catalog as
;; an association list
;    
;    if keyword_set(gr_hstdetect) then begin
;
;       assocname_hst = catpath+'sg1120_hst_merged_gr_07apr27.assoc'
;       if (n_elements(hstcat) eq 0L) then begin
;          catname_hst = catpath+'sg1120_hst_merged_07apr27.cat'
;          splog, 'Reading '+catname_hst
;          hstcat = rsex(catname_hst)
;          hdr = headfits(mosaicpath+'sg1120_rprime.fits') ; g and r have the same astrometry header
;          extast, hdr, astr
;          ad2xy, hstcat.alpha_j2000, hstcat.delta_j2000, astr, hst_gr_x, hst_gr_y
;          splog, 'Writing '+assocname_hst
;          openw, lun, assocname_hst, /get_lun
;          for ii = 0L, n_elements(hstcat)-1L do printf, lun, ii, hst_gr_x[ii], $
;            hst_gr_y[ii], format='(I6,2F12.4)'
;          free_lun, lun
;       endif
;
;       pixscale = 0.188
;       assoc_radius = strtrim(0.5/pixscale,2) ; 0.5 arcsec [pixel]
;       
;       photaper = [1.0,2.0,3.0] ; [diameter,arcsec]
;       photaper = strjoin(strtrim(string(photaper/pixscale,format='(F10.2)'),2),',') ; [diameter,pixel]
;    
;       seeing = '0.9'  ; figure this out iteratively
;       magzero = '0.0' ; string(sxpar(hdr,'MAGZERO'),format='(G0.0)')
;
;; detect on the r-band image, measure on the r-band image       
;       
;       imagelist = mosaicpath+['sg1120_rprime.fits','sg1120_rprime.fits']
;       weightlist = repstr(imagelist,'.fits','.weight.fits')
;       catname = catpath+'sg1120_rprime_hstmatched_'+date+'.cat'
;
;       checkimage_type = 'BACKGROUND,SEGMENTATION,APERTURES'
;       checkimage_name = mosaicpath+strjoin(['rprime.background.fits','rprime.segmentation.fits','rprime.apertures.fits'],',')
;
;       spawn, 'sex '+strjoin(imagelist,',')+' -c '+sexconfig+' -PARAMETERS_NAME '+sexparam+$
;         ' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -CATALOG_NAME '+catname+' -CATALOG_TYPE ASCII_HEAD'+$
;         ' -DETECT_THRESH 1.5 -DETECT_MINAREA 5 -DEBLEND_MINCONT 0.003'+$
;         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+strjoin(weightlist,',')+' -WEIGHT_THRESH 0 -WEIGHT_GAIN N'+$
;         ' -PHOT_APERTURES '+photaper+' -MAG_ZEROPOINT '+magzero+' -SEEING_FWHM '+seeing+$
;         ' -ASSOC_NAME '+assocname_hst+' -ASSOC_PARAMS 2,3 -ASSOC_DATA 1,2,3 -ASSOC_RADIUS '+assoc_radius+' -ASSOCSELEC_TYPE MATCHED -ASSOC_TYPE NEAREST'+$
;         ' -MEMORY_OBJSTACK 30000 -MEMORY_PIXSTACK 3000000 -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 2 '+$
;         ' -CHECKIMAGE_TYPE '+checkimage_type+' -CHECKIMAGE_NAME '+checkimage_name, /sh
;
;; detect on the r-band image, measure on the g-band image       
;       
;       imagelist = mosaicpath+['sg1120_rprime.fits','sg1120_gprime.fits']
;       weightlist = repstr(imagelist,'.fits','.weight.fits')
;       catname = catpath+'sg1120_gprime_hstmatched_'+date+'.cat'
;
;       checkimage_type = 'BACKGROUND,SEGMENTATION,APERTURES'
;       checkimage_name = mosaicpath+strjoin(['gprime.background.fits','gprime.segmentation.fits','gprime.apertures.fits'],',')
;
;       spawn, 'sex '+strjoin(imagelist,',')+' -c '+sexconfig+' -PARAMETERS_NAME '+sexparam+$
;         ' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+' -CATALOG_NAME '+catname+' -CATALOG_TYPE ASCII_HEAD'+$
;         ' -DETECT_THRESH 1.5 -DETECT_MINAREA 5 -DEBLEND_MINCONT 0.003'+$
;         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+strjoin(weightlist,',')+' -WEIGHT_THRESH 0 -WEIGHT_GAIN N'+$
;         ' -PHOT_APERTURES '+photaper+' -MAG_ZEROPOINT '+magzero+' -SEEING_FWHM '+seeing+$
;         ' -ASSOC_NAME '+assocname_hst+' -ASSOC_PARAMS 2,3 -ASSOC_DATA 1,2,3 -ASSOC_RADIUS '+assoc_radius+' -ASSOCSELEC_TYPE MATCHED -ASSOC_TYPE NEAREST'+$
;         ' -MEMORY_OBJSTACK 30000 -MEMORY_PIXSTACK 3000000 -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 2 '+$
;         ' -CHECKIMAGE_TYPE '+checkimage_type+' -CHECKIMAGE_NAME '+checkimage_name, /sh
;
;    endif       
       

return
end
    
