pro ediscs_specfit, index=index1, solar=solar, test=test, doplot=doplot
; jm06sep27nyu
; jm08may20nyu - 
; jm09may22nyu - major rewrite
; jm09aug18ucsd - v2.1; solar metallicity is the default; also fit
;   H-delta in emission now
; jm09oct02ucsd - do not cross-correlate

    zsnrcross = -1.0 ; do not cross-correlate
    solar = 1 ; ASSUME SOLAR METALLICITY
    
    if keyword_set(doplot) or keyword_set(test) then $
      silent = 0 else silent = 1

    spec1dpath = ediscs_path(/spec1d)
    version = ediscs_version(/specfit)

    base_specfitpath = ediscs_path(/specfit)
    specfitpath = base_specfitpath+version+'/'

    linefile = base_specfitpath+'elinelist_'+version+'.dat'
    indexfile = base_specfitpath+'indexlist_'+version+'.dat'
    if keyword_set(solar) then Zbc03 = ['Z02'] else $
      Zbc03 = ['Z004','Z02','Z05']

    templatefile = base_specfitpath+'BC03_'+Zbc03+$
      '_chabrier_ediscs_templates.fits'

; fit each cluster separately
    ediscs = read_ediscs(/ancillary)

    allcluster = strtrim(ediscs.cluster,2)
    cluster = allcluster[uniq(allcluster,sort(allcluster))]
    ncluster = n_elements(cluster)

    specres = 5.5
    
    t0 = systime(1)
;   for iclust = 0, 0 do begin
    for iclust = 0, ncluster-1 do begin
;   for iclust = ncluster/2L, ncluster-1L do begin

       splog, 'Fitting cluster '+cluster[iclust]
       suffix = 'ediscs_'+cluster[iclust]
       if keyword_set(test) then suffix = suffix+'_test'

       these = where(cluster[iclust] eq allcluster,nobj)
       if (n_elements(index1) ne 0) then index = index1 else $
         index = lindgen(nobj)
       nindex = n_elements(index)

       npix = 1E4 & for jj = 0, nindex-1 do npix = npix < sxpar(headfits($
         spec1dpath+strtrim(ediscs[these[index[jj]]].specfile,2)),'NAXIS1')
       wave = fltarr(npix,nindex)
       flux = fltarr(npix,nindex)
       invvar = fltarr(npix,nindex)
       for jj = 0, nindex-1 do begin
          thisfile = strtrim(spec1dpath+ediscs[these[index[jj]]].specfile,2)
          hh = headfits(thisfile)
          npix1 = sxpar(hh,'NAXIS1')
          flux[*,jj] = (mrdfits(thisfile,0,/silent))[0:npix-1]
          invvar[*,jj] = (1.0/mrdfits(thisfile,1,/silent)^2)[0:npix-1]
          wave[*,jj] = (make_wave(hh))[0:npix-1]
       endfor

       zobj = ediscs[these[index]].z
       galaxy = ediscs[these[index]].galaxy

       specdata = ispeclinefit(wave,flux,invvar,zobj=zobj,specres=specres,$
         vdisp=vdisp,sigmax=sigmax,zsnrcross=zsnrcross,galaxy=galaxy,$
         outpath=specfitpath,suffix=suffix,templatefile=templatefile,$
         linefile=linefile,indexfile=indexfile,specfit=specfit,/nologfile,$
         vmaxtweak=vmaxtweak,vlinemaxtweak=vlinemaxtweak,/clobber,$
         silent=silent,doplot=doplot)
             
;      specdata = ispeclinefit(speclist,specres=5.9,snrcut=0.0,dustmodel=0,zsnrcross=1.0,$
;        datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=0,sigmax=800.0,$
;        /charlot,/zcrosscor,postscript=postscript,/write,vmaxshift=400.0,nologfile=1,$
;        starvdisp=100.0,doplot=doplot,eigenfile=eigenfile,eigendir=eigendir);,/positive_emission)

    endfor
    splog, 'Total '+string((systime(1)-t0)/3600.0,format='(G0.0)')+' hours'

    ediscs_merge_specfit, /write

return
end    
