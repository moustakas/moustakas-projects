pro hst_merge_stefan_catalogs, bigcat, finalcat, match1=match1, match2=match2, dcr=dcr, write=write
; jm07apr26nyu - merge the catalogs of the 10 independent HST
;                pointings generated by Stefan
; jm07jun25nyu - merge the catalogs generated off the Gonzalez
;                reductions 

    date = strjoin(strsplit(strlowcase(strmid(im_today(),2,9)),' ',/ext))
;   get_date, date, timetag=0 & date = repstr(date,'-','_')

    if (n_elements(dcr) eq 0L) then dcr = 0.5D ; [arcsec]
    
    datapath = hst_path()+'gonzalez_drizzle/'
    catpath = hst_path(/catalogs)+'gonzalez_deblend_0.015/'
;   catpath = hst_path()+'catalogs/'
;   catpath = sg1120_path(/catalogs)
;   catpath = hst_path(/catalogs)+'blakeslee_deblend_0.009/'

    if (n_elements(bigcat) eq 0L) then begin
       imagelist = file_search(datapath+'hst_f814w_drz_pos??.fits')
       catlist = [file_search(catpath+'SG_p?.cat'),catpath+'SG_p10.cat']
;      catlist = [file_search(catpath+'sexcat_pos?.cat'),catpath+'sexcat_pos10.cat']
       xyreglist = repstr(catlist,'.cat','.xy.reg')
       radecreglist = repstr(catlist,'.cat','.radec.reg')
       ncat = n_elements(catlist)
       for icat = 0L, ncat-1L do begin
          print, format='("Reading catalog ",I0,"/",I0,".",A10,$)', icat+1, ncat, string(13b)
          thiscat = rsex(catlist[icat],use_row=1L) ; read the first row because FLUX_APER is LONG, not DOUBLE
          nobj = n_elements(thiscat)
; ---------------------------------------------------------------------------
; NOTE!  For some reason, the _WORLD coordinates are wrong/off.
; Possibly a bug in SE?
          hdr = headfits(imagelist[icat],ext=1L)
          extast, hdr, astr
          xy2ad, thiscat.xwin_image, thiscat.ywin_image, astr, aa, dd
          thiscat.xwin_world = aa & thiscat.ywin_world = dd
; ---------------------------------------------------------------------------
          struct_print, struct_trimtags(thiscat,select=['XWIN_IMAGE','YWIN_IMAGE'],$
            format=['G15.8','G15.8']), file=xyreglist[icat], /no_head
          struct_print, struct_trimtags(thiscat,select=['XWIN_WORLD','YWIN_WORLD'],$
            format=['E16.10','E17.10']), file=radecreglist[icat], /no_head
          i1 = strpos(catlist[icat],'p')+1 & i2 = strpos(catlist[icat],'.cat')
;         i1 = strpos(catlist[icat],'pos')+3 & i2 = strpos(catlist[icat],'.cat')
          pointing = fix(strmid(catlist[icat],i1,i2-i1)); & print, pointing
          moretags = replicate({pointing: fix(pointing), original_number: 0L, nmulti: 0L},nobj)
          cat1 = struct_addtags(thiscat,moretags)
          cat1.original_number = thiscat.number
          if (icat eq 0L) then bigcat = cat1 else bigcat = [bigcat,cat1]
          help, cat1
       endfor
    endif
    ngalaxy = n_elements(bigcat)

; identify duplicates; among the duplicates go through and keep
; objects that are on the same pointing (NMULTI_SAME), since these are
; likely to be deblended from the same physical object (e.g., HII
; regions, etc.); for duplicates (multiples) coming from different
; pointings keep the one with the highest S/N magnitude measurement 

    splog, 'Running friends-of-friends using DCR='+string(dcr,format='(F3.1)')+' arcsec.'
    ingroup = spheregroup(bigcat.xwin_world,bigcat.ywin_world,dcr/3600D0,$
      multgroup=mult,firstgroup=first,nextgroup=next)
    ugroup = ingroup[uniq(ingroup,sort(ingroup))] ; unique groups

    finalindx = lindgen(ngalaxy)
    for ii = 0L, n_elements(ugroup)-1L do begin
       match = where((ingroup eq ugroup[ii]),nmatch)
       point = bigcat[match].pointing
       if (nmatch gt 1) then begin
;         struct_print, /no_head, struct_trimtags(bigcat[match],select=['POINTING',$
;           'ORIGINAL_NUMBER','XWIN_WORLD','YWIN_WORLD','MAG_AUTO','MAGERR_AUTO','FLUX_RADIUS'])
          upoint = point[uniq(point,sort(point))]
          nupoint = n_elements(upoint)
          case nmatch of
; cases: (1) same pointing (over-deblended source; pick the
; brightest); (2) different pointing (true duplicate: select on S/N) 
             2L: begin
                case nupoint of
                   1L: begin
                      bigcat[match].nmulti = nmatch ; case (1) [potentially dumb]
                      minmag = min(bigcat[match].mag_auto,minmagindx)
                      remove, minmagindx, match & finalindx[match] = -1L
                   end
                   2L: begin
                      bigcat[match].nmulti = nmatch ; case (2)
                      maxsnr = max(bigcat[match].flux_best/bigcat[match].fluxerr_best,maxsnrindx)
;                     maxsnr = max(bigcat[match].mag_auto/bigcat[match].magerr_auto,maxsnrindx)
                      remove, maxsnrindx, match & finalindx[match] = -1L
                   end
                endcase
             end
; cases: (1) all from the same pointing (over-deblending source)
; (although in some cases the magnitudes are all comparable!!; pick
; the brightest one); (2) two objects from the same pointing and one
; from another pointing (pick the brightest); (3) all three from a
; different pointing (true multiple: select on S/N)
             3L: begin
                case nupoint of
                   1L: begin
                      bigcat[match].nmulti = nmatch ; case (1) [potentially dumb]
                      minmag = min(bigcat[match].mag_auto,minmagindx) 
                      remove, minmagindx, match & finalindx[match] = -1L
                   end
                   2L: begin
                      bigcat[match].nmulti = nmatch ; case (2)
                      minmag = min(bigcat[match].mag_auto,minmagindx)
                      remove, minmagindx, match & finalindx[match] = -1L
                   end
                   3L: begin
                      bigcat[match].nmulti = nmatch ; case (3)
                      maxsnr = max(bigcat[match].flux_best/bigcat[match].fluxerr_best,maxsnrindx)
;                     maxsnr = max(bigcat[match].mag_auto/bigcat[match].magerr_auto,maxsnrindx)
                      remove, maxsnrindx, match & finalindx[match] = -1L
                   end
                endcase
             end
; case: (1) *pairs* of objects from different pointings (pick the
; first pair, and then pick the brightest)
             4L: begin
                bigcat[match].nmulti = nmatch ; case (3)
                toss = where(point ne upoint[0],comp=keep) ; keep the first pair!
                finalindx[match[toss]] = -1L
                minmag = min(bigcat[match[keep]].mag_auto,minmagindx) ; pick the brightest [potentially dumb]
                toss2 = match[keep] & remove, minmagindx, toss2 & finalindx[toss2] = -1L
             end
             else: begin ; keep the brightest one
                bigcat[match].nmulti = nmatch
                maxmag = max(bigcat[match].mag_auto,maxmagindx)
                remove, maxmagindx, match & finalindx[match] = -1L
                finalindx[match] = -1L ; throw out?!
             end
          endcase
       endif 
    endfor

    mult = where((finalindx eq -1L),nmult,comp=keep,ncomp=nkeep)
    finalcat = bigcat[finalindx[keep]]
    finalcat.number = lindgen(nkeep)+1L
    
    splog, 'Identified and removed '+string(nmult,format='(I0)')+' multiples'

; test code:

;    splog, 'Running test code.'
;    ingroup = spheregroup(finalcat.xwin_world,finalcat.ywin_world,dcr/3600D0)
;    ugroup = ingroup[uniq(ingroup,sort(ingroup))] ; unique groups
;
;    if (nkeep ne n_elements(ugroup)) then begin
;       for ii = 0L, n_elements(ugroup)-1L do begin
;          match = where((ingroup eq ugroup[ii]),nmatch)
;          if (nmatch gt 1) then begin
;             struct_print, /no_head, struct_trimtags(finalcat[match],select=['POINTING',$
;               'ORIGINAL_NUMBER','XWIN_WORLD','YWIN_WORLD','MAG_AUTO','MAGERR_AUTO','FLUX_RADIUS'])
;          endif
;       endfor
;    endif
    
    if keyword_set(write) then begin
       catfile = catpath+'sg1120_hst.cat'
;      catfile = catpath+'sg1120_hst_dcr'+string(dcr,format='(F3.1)')+'.cat'
       splog, 'Writing '+catfile
;      w = where(bigcat.pointing eq 10L)
;      struct_print, struct_trimtags(bigcat[w],select=['XWIN_WORLD','YWIN_WORLD'],$
;        format=['E16.10','E17.10']), file=repstr(catfile,'.cat','.radec.reg'), /no_head
       struct_print, struct_trimtags(finalcat,select=['XWIN_WORLD','YWIN_WORLD'],$
         format=['E16.10','E17.10']), file=repstr(catfile,'.cat','.radec.reg'), /no_head
       wsex, finalcat, outfile=catfile
    endif
    
return
end

;;; cross-match the catalog with itself to identify duplicates    
;;
;;    if (n_elements(match1) eq 0L) and (n_elements(match2) eq 0L) then begin
;;       splog, 'Matching catalogs.'
;;       spherematch, bigcat.xwin_world, bigcat.ywin_world, bigcat.xwin_world, $
;;         bigcat.ywin_world, dcr/3600D0, match1, match2, dist12, maxmatch=1
;;;      spherematch, bigcat.alpha_j2000, bigcat.delta_j2000, bigcat.alpha_j2000, $
;;;        bigcat.delta_j2000, dcr/3600.0, match1, match2, dist12, maxmatch=1
;;    endif
;;
;;    splog, 'Removing duplicates.'
;;    hist = histogram(match1,binsize=1,reverse_indices=indx)
;;    bigcat.nmulti = fix(hist-1) ; exclude the object itself
;;    finalindx = lonarr(ngalaxy)
;;
;;    for igal = 0L, ngalaxy-1L do begin
;;;      print, format='("Galaxy ",I0,"/",I0,".",A10,$)', igal+1, ngalaxy, string(13b)
;;       thisindx = match2[indx[indx[igal]:indx[igal+1]-1]] ; this has to be MATCH2
;;       if (hist[igal] eq 1L) then finalindx[igal] = thisindx else begin
;;; require "duplicates" to belong to different pointings
;;          thispointing = bigcat[thisindx].pointing
;;;          if n_elements(thisindx) gt 2L then stop
;;          if monotonic(thispointing[sort(thispointing)]) then begin
;;             dup_info1 = struct_trimtags(bigcat[thisindx],select=['POINTING','ORIGINAL_NUMBER',$
;;               'XWIN_WORLD','YWIN_WORLD','MAG_AUTO','MAGERR_AUTO'])
;;             if (n_elements(dup_info) eq 0L) then dup_info = dup_info1 else $
;;               dup_info = [dup_info,dup_info1]
;;             magsnr = bigcat[thisindx].mag_auto/bigcat[thisindx].magerr_auto
;;             maxsnr = max(magsnr,maxsnrindx) ; select by S/N
;;             finalindx[igal] = thisindx[maxsnrindx]
;;          endif
;;       endelse
;;    endfor
;;
;;    finalindx = finalindx[rem_dup(finalindx)] ; remove duplicates
;;    finalcat = bigcat[finalindx]
;;    finalcat.number = lindgen(n_elements(finalcat))+1L

