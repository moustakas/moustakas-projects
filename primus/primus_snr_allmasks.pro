function primus_snr_summary, date, basedir
; read all the SNR files
    
    for ii = 0L, n_elements(date)-1L do begin
       thisdir = basedir+date[ii]+'/'
       qaplotsdir = thisdir+'qaplots/'
       snrdir = qaplotsdir+'snr/'
       allmask = file_basename(file_search(thisdir+'*_summary.txt',count=nmask))
       for jj = 0L, nmask-1L do begin
          thismask = strmid(allmask[jj],0,strpos(allmask[jj],'_summary.txt'))
; read the summary file to figure out the exposures
          info = djs_readlines(thisdir+allmask[jj])
          info = info[where(strmatch(info,'*Extraction Statistics*'))+1:$
            where(strmatch(info,'*Signal/Noise Statistics*'))-1]
          info = info[where(strcompress(info,/remove) ne '')]
;         niceprint, kinfo
          openw, lun, '/tmp/junk', /get_lun
          for mm = 0L, n_elements(info)-1L do $
            printf, lun, info[mm]
          free_lun, lun
          readcol, '/tmp/junk', ccd, skipline=1, format='A', /silent
          nccd = n_elements(ccd)
;         niceprint, ccd
          ccd_start = repstr(ccd[0],'ccd','')
          ccd_end = repstr(ccd[nccd-1],'ccd','')
; build the list of parameter files
          snrfile_mask = snrdir+thismask+'_snrfits.par'
          snrfile_ccd = snrdir+ccd+'_snrfits.par'
          nsnr = n_elements(snrfiles)
;         niceprint, snrfiles
;         if (total(file_test(snrfiles,/regular)) ne nsnr) then message, 'Expected files missing'
; read all the parameter files and store the results
          snr_mask = yanny_readone(snrfile_mask)
          summary1 = {date: date[ii], mask: thismask, snr_mask: 0.0, ccd: '', snr_ccd: 0.0}
          summary1 = replicate(summary1,nccd)
          summary1.snr_mask = 10.0^poly(0.0,snr_mask.magfit)
          for kk = 0L, nccd-1L do begin
             snr_ccd = yanny_readone(snrfile_ccd[kk])
             summary1[kk].ccd = ccd[kk]
             summary1[kk].snr_ccd = 10.0^poly(0.0,snr_ccd.magfit)
          endfor
;         struct_print, summary1
          if (n_elements(summary) eq 0L) then summary = summary1 else $
            summary = [summary,summary1]
       endfor
    endfor
    
return, summary
end

pro primus_snr_allmasks, date=date1
; jm09jan06nyu - plot the S/N of all the PRIMUS masks

    rerun = '0020'
    basedir = getenv('PRIMUS_DATA')+'/redux/'+rerun+'/'

    summary_allfile = basedir+'snr_summary_'+rerun+'_all.fits'
    summary_maskfile = basedir+'snr_summary_'+rerun+'_masks.fits'
    summary_qaplot = basedir+'snr_summary_'+rerun+'.ps'
    
    if (n_elements(date1) eq 0L) then $
      date = file_basename(file_search(basedir+'ut??????')) else date = date1

; read all the SNR files    
    sum = primus_snr_summary(date,basedir)

    allmask = strtrim(sum.mask,2)
    uu = uniq(allmask,sort(allmask))
    srt = sort(sum[uu].snr_mask)
;   srt = sort(sum[uu].date)
    usum = sum[uu[srt]]
    nusum = n_elements(usum)
;   struct_print, usum

    umask = allmask[uu]
    numask = n_elements(umask)

; write out    
    splog, 'Writing '+summary_allfile
    mwrfits, sum, summary_allfile, /create
    splog, 'Writing '+summary_maskfile
    mwrfits, usum, summary_maskfile, /create

    splog, 'Writing '+summary_qaplot
    im_plotfaves, /post
    dfpsplot, summary_qaplot, /color, /landscape
    datecolors = ['red','blue','dark green']
; page 1
    djs_plot, [0], [0], /nodata, /xsty, ysty=3, xrange=[0,numask+1], $
      yrange=[min(sum.snr_mask)<min(sum.snr_ccd),max(sum.snr_mask)>max(sum.snr_ccd)], $
      xtitle='Mask Number', ytitle='Fiducial S/N'
    for jj = 0L, numask-1L do begin
       these = where((umask[jj] eq allmask),nthese)
; individual exposures
       if (nthese eq 1L) then plots, replicate(jj+1L,nthese), sum[these].snr_ccd, psym=7 else $
         djs_oplot, replicate(jj+1L,nthese), sum[these].snr_ccd, psym=7, sym=1.5
; coadded masks
       alldate = sum[these].date
       udate = alldate[uniq(alldate,sort(alldate))]
       for kk = 0L, n_elements(udate)-1L do begin
          ww = where(udate[kk] eq alldate)
          plots, jj+1, sum[these[ww[0]]].snr_mask, psym=4, sym=3, $
            color=djs_icolor(datecolors[kk]), thick=5.0
       endfor
    endfor
;;; page 2
;;    djs_plot, [0], [0], /nodata, /xsty, ysty=3, xrange=[0,nusum+1], $
;;      yrange=minmax(usum.snr_mask), xtitle='UT Date', ytitle='Fiducial S/N';, $
;;;     xtickname=[' ',repstr(usum.date,'ut','')], xtickinterval=2, xminor=-1
;;;   djs_oplot, lindgen(nusum)+1L, usum.snr_mask, psym=-4, sym=3, color='red'
;;; coadded masks
;;    alldate = usum.date
;;    udate = alldate[uniq(alldate,sort(alldate))]
;;    datecounter = 1
;;    for kk = 0L, n_elements(udate)-1L do begin
;;       ww = where(udate[kk] eq alldate,nww)
;;stop
;;       if (nww gt 1L) then oplot, replicate(datecounter,nww), usum[jj[ww]].snr_mask, $
;;         psym=4, sym=3, color=djs_icolor(datecolors[kk]), thick=5.0 else $
;;       plots, datecounter, usum[jj[ww[0]]].snr_mask, psym=4, sym=3, $
;;         color=djs_icolor(datecolors[kk]), thick=5.0
;;       datecounter = datecounter + 1L
;;       stop
;;    endfor
    dfpsclose
    im_plotfaves
    
stop    
    
return
end
