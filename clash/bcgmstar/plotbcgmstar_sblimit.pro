pro plotbcgmstar_sblimit
; jm14may15siena - plot the SB limits as a function of cluster and
; filter wavelength 

    paperpath = bcgmstar_path(/paper)
    skyinfopath = bcgmstar_path()+'skyinfo/'

    sample = read_bcgmstar_sample()
;   sample = sample[sort(sample.mvir)]
    cluster = strtrim(sample.shortname,2)
    ncl = n_elements(sample)

    filt = bcgmstar_filterlist(short=short,weff=weff,pivot=pivot)
    nfilt = n_elements(filt)
    
; gather the data for the box-and-whisker plot
    data = fltarr(nfilt,ncl)-1.0
    for ic = 0, ncl-1 do begin
       sky = mrdfits(skyinfopath+'skyinfo-'+cluster[ic]+'.fits.gz',1,/silent)
       for jj = 0, n_elements(sky)-1 do begin
          this = where(strtrim(sky[jj].band,2) eq strtrim(short,2))
          data[this,ic] = sky[jj].sblimit-2.5*alog10(3.0) ; 3-sigma
       endfor
    endfor
    num = fix(total(data gt 0,2))
    niceprint, num, filt
    
    psfile = paperpath+'bcgmstar_sblimit.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, charsize=1.6, ymargin=[1.0,1.1]

    cgboxplot, data, xtitle='', $ ; textoidl('Filter Wavelength (\mu'+'m)'), $
      ytitle=textoidl('Surface Brightness Limit (3\sigma, mag arcsec^{-2})'), $
      /fillboxes, boxcolor='plum', AxisColor='black', stats=stats, $
      OutLineColor='navy', OutlierColor='grn5', $ ; xlocation=weff/1D4, $
      missing_data=-1.0, position=pos, yrange=[25.5,22], $ ; xrange=[0.3,1.6], 
      labels=strupcase(short), rotate=-45, width=0.6, ytickinterval=1
;     output=paperpath+'bcgmstar_sblimit.pdf'
    im_plotconfig, psfile=psfile, /psclose, /pdf

;   struct_print, stats
;   print, weighted_quantile(data[0,*],quant=[0.25,0.5,0.75])
    niceprint, short, weff, pivot, stats.median

stop    
    
return
end
    

