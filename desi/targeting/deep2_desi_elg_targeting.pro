function get_dndm, rmag, weight=weight, faintcut=faintcut, $
  brightcut=brightcut, magaxis=magaxis
; get the number counts    
    area = 0.4342 ; deg^2
    binsize = 0.2
    nbins = ceil((faintcut-brightcut)/binsize)
    magaxis = lindgen(nbins)*binsize+brightcut+binsize/2.0
    dndm = hogg_histogram(rmag,[brightcut,faintcut],$ ; #/0.2 mag/deg^2
      nbins,weight=weight)/area
;   dndm = im_hist1d(rmag,weight,min=brightcut,$ 
;     max=faintcut,binsize=binsize)/area/binsize
return, dndm
end

function get_hiz, ugriz, mostek=mostek
; select high-redshift galaxies using a color-cut
    gr = ugriz[1,*]-ugriz[2,*]
    rz = ugriz[2,*]-ugriz[4,*]
; fiducial cut    
    rzcut = 1.2
    blue = where(rz lt rzcut,comp=red,nblue,ncomp=nred)
    if nblue ne 0 then hiz_blue = where(gr[blue] lt (0.7*rz[blue]+0.05))
    if nred ne 0 then hiz_red = where(gr[red] lt (1.4*rz[red]-0.79))
    if nblue ne 0 and nred ne 0 then hiz = [blue[hiz_blue],red[hiz_red]]
    if nblue ne 0 and nred eq 0 then hiz = blue[hiz_blue]
    if nblue eq 0 and nred ne 0 then hiz = red[hiz_red]
; Mostek's cut    
    if keyword_set(mostek) then hiz = where(gr lt (0.68*rz-0.08) and $
      (rz gt 0.2) and (rz lt 1.3))
return, hiz
end

function get_ugriz, cat, unwise=unwise
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii, unwise=unwise
    cat = struct_addtags(cat,replicate({ugriz: fltarr(5), $
      ugriz_err: fltarr(5), bri: fltarr(3), bri_err: fltarr(3), $
      wise: fltarr(2), wise_err: fltarr(2)},$
      n_elements(cat)))
    cat = im_struct_trimtags(cat,select=tag_names(cat),newtags=$
      repstr(repstr(tag_names(cat),'RA_DEEP','ra'),'DEC_DEEP','DEC'))
    cat = struct_trimtags(cat,except=['RA_SDSS','DEC_SDSS'])
    mag = maggies2mag(mm,ivarmaggies=ii,magerr=magerr)
    cat.bri = mag[0:2,*]
    cat.bri_err = magerr[0:2,*]
    cat.ugriz = mag[3:7,*]
    cat.ugriz_err = magerr[3:7,*]
    if keyword_set(unwise) then begin
       cat.wise = mag[8:9,*]
       cat.wise_err = magerr[8:9,*]
    endif
return, cat
end

pro deep2_desi_elg_targeting, build_parent=build_parent, qaplots_parent=qaplots_parent, $
  build_oiihiz=build_oiihiz, qaplots_select=qaplots_select
; jm14mar01siena - optimize the DESI targeting using DEEP2
; spectroscopy of galaxies in the Groth Strip (Field 1)

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/v1.1/'
    isedfit_paramfile = templatepath+'desi_deep2_paramfile.par'
    winpath = deep2_path(/window)
    outpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'
    technotepath = getenv('IM_SVNREPOS')+'/desi/technotes/targeting/elg-deep2/trunk/'

    area = 0.4342 ; deg^2
;   area = 0.6 ; [deg^2]

    brightcut = 18.5
    faintcut = 24.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    
; ###########################################################################
; build a parent photometric sample concentrated just in the EGS with
; some basic sample cuts 
    if keyword_set(build_parent) then begin
       allphot = mrdfits(deep2_path(/cat)+'deep2.pcat_ext.fits.gz',1)    

       egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
       phot_egs = get_ugriz(allphot[egs])
       pointing_egs = strmid(strtrim(phot_egs.objno,2),0,2) ; test

; select objects in the spectroscopic footprint in a reasonable
; magnitude range with PGAL>0.5
       win = mrdfits(winpath+'windowf.egs.fits.gz',0,hdr)
       extast, hdr, astr
       area = total(win gt 0)*determ(astr.cd) 
       splog, 'Spectroscopic area (deg^2) = ', area
       
       ad2xy, phot_egs.ra, phot_egs.dec, astr, xx, yy
       maskweight = interpolate(win,xx,yy,missing=0)

       keep = where(maskweight gt 0.25 and phot_egs.pgal ge 0.5 and $
         phot_egs.badflag eq 0 and $
         (phot_egs.bestr + 2.5*alog10(!pi*(3*phot_egs.rg*0.207)^2)) le 26.5 and $
         phot_egs.ugriz[2,*] gt brightcut and phot_egs.ugriz[2,*] lt faintcut,nphot)
       phot = phot_egs[keep]

       nmissing = total(phot.ugriz[2,*] eq 0) ; missing r-band photometry --> none!
       splog, 'Parent photometric sample = ', nphot, negs
       splog, 'Missing r-band photometry = ', nmissing

       im_mwrfits, phot, outpath+'egs_parent_phot.fits', /clobber

; also write out a sample of stars so we can investigate the stellar
; locus
       star = where(phot_egs.pgal lt 0.2 and phot_egs.badflag eq 0 and $
         phot_egs.ugriz[2,*] gt brightcut and phot_egs.ugriz[2,*] lt faintcut,nstar)
       splog, 'Sample of stars', nstar
       im_mwrfits, phot_egs[star], outpath+'egs_stars_phot.fits', /clobber
       
; now build our parent spectroscopic sample of galaxies in the EGS
; field (apply the cut on Q>=3 below)
       allzcat = read_deep2_zcat(phot=allzcat_phot,weight=allzcat_weight,/all)
       allzcat = struct_addtags(struct_addtags(allzcat,struct_trimtags(allzcat_weight,select=$
         ['obj_weight','targ_weight'])),replicate({mask_weight: -1.0},n_elements(allzcat)))
       egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',negs)

       zcat_egs = allzcat[egs]
       zcat_phot_egs = get_ugriz(allzcat_phot[egs],/unwise)

; this cut is to choose objects that are within the spectroscopic
; footprint of Field 1       
       pointing1 = strmid(strtrim(zcat_egs.objno,2),0,2)
       ad2xy, zcat_egs.ra, zcat_egs.dec, astr, xx, yy
       zcat_maskweight1 = interpolate(win,xx,yy,missing=0)
       zcat_egs.mask_weight = zcat_maskweight1

       keep = where(zcat_maskweight1 gt 0 and zcat_egs.targ_weight gt 0 and $
         zcat_phot_egs.ugriz[2,*] gt brightcut and zcat_phot_egs.ugriz[2,*] lt faintcut,nspec)
       zcat = zcat_egs[keep]
       zcat_phot = zcat_phot_egs[keep]

       nspec_weighted = long(total(1.0/zcat.targ_weight))
       nmissing = total(zcat_phot.ugriz[2,*] eq 0) ; missing r-band photometry --> none!
       splog, 'Parent spectroscopic sample = ', negs, nspec, nspec_weighted
       splog, 'Missing r-band photometry = ', nmissing

       im_mwrfits, zcat, outpath+'egs_parent_zcat.fits', /clobber
       
; also write out a Q>=3 sample with the corresponding [OII] fluxes and
; upper limits 
       ppxf = read_deep2(/ppxf,/fixoii)
       match, zcat.objno, ppxf.objno, m1, m2
       srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
       zcat_phot_q34 = zcat_phot[m1]
       zcat_q34 = zcat[m1]
       ppxf_q34 = ppxf[m2]

; build the [OII] flux       
       kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
         'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1,$
         rows=m2)
       oii = deep2_get_oiiflux(ppxf_q34,cflux_3727_rest=kised.cflux_3727)
       zcat_q34 = struct_addtags(struct_addtags(zcat_q34,struct_trimtags(zcat_phot_q34,$
         select=['ugriz','ugriz_err','bri','bri_err','wise','wise_err'])),oii)
       
       im_mwrfits, zcat_q34, outpath+'egs_parent_zcat.Q34.fits', /clobber
    endif

; ###########################################################################
; read the parent sample and build a sample of high-redshift, strong
; [OII]-emitting galaxies     
    if keyword_set(build_oiihiz) then begin
       zcat = mrdfits(outpath+'egs_parent_zcat.fits.gz',1)
       kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
         'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1,$
         rows=zcat.parentindx)

       ppxf = (read_deep2(/ppxf))[zcat.parentindx]
       ngal = n_elements(zcat)

       colorhiz = get_hiz(zcat.ugriz)
       colorhiz_mostek = get_hiz(zcat.ugriz,/mostek)
       ncolorhiz = n_elements(colorhiz)
       
       hiz = where(zcat.z gt 0.7,nhiz,comp=loz,ncomp=nloz)
       oii = deep2_get_oiiflux(ppxf,kised,oiilimit=oiilimit)

; assume that the objects with formal flux limits above our [OII] cut
; are detections (this is a small number of objects)    
       oiibright = where(oii[colorhiz] ge oiicut1)
       oiifaint = where(oii[colorhiz] lt oiicut1 and oii[colorhiz] ne -2.0)
       oiinone = where(oii[colorhiz] eq -2.0)

       dndm = get_dndm(zcat.ugriz[2,*],weight=1.0/zcat.weight,$
         faintcut=faintcut,brightcut=brightcut,magaxis=magaxis)
       dndm_colorhiz = get_dndm(zcat[colorhiz].ugriz[2,*],$
         weight=1.0/zcat[colorhiz].weight,$
         faintcut=faintcut,brightcut=brightcut)
       dndm_oiibright = get_dndm(zcat[colorhiz[oiibright]].ugriz[2,*],$
         weight=1.0/zcat[colorhiz[oiibright]].weight,$
         faintcut=faintcut,brightcut=brightcut)

; get the *ratio* of the number of bright-to-faint [OII] sources so
; that we can correct for the objects with missing [OII] measurements
; (see plot below)     
       dndm_oiifaint = get_dndm(zcat[colorhiz[oiifaint]].ugriz[2,*],$
         weight=1.0/zcat[colorhiz[oiifaint]].weight,$
         faintcut=faintcut,brightcut=brightcut)
       dndm_oiinone = get_dndm(zcat[colorhiz[oiinone]].ugriz[2,*],$
         weight=1.0/zcat[colorhiz[oiinone]].weight,$
         faintcut=faintcut,brightcut=brightcut)
       denom = dndm_oiibright+dndm_oiifaint
       dndm_oiicor = dndm_oiinone*dndm_oiibright/(denom+(denom eq 0))*(denom ne 0)
       dndm_oiibright_cor = dndm_oiibright + dndm_oiicor    
    endif
       
; ###########################################################################
; build QAplots of our high-z target selection
    if keyword_set(qaplots_select) then begin
;      phot = mrdfits(outpath+'egs_parent_phot.fits.gz',1)
;      zcat = mrdfits(outpath+'egs_parent_zcat.fits.gz',1)
       zcat = mrdfits(outpath+'egs_parent_zcat.Q34.fits.gz',1)

       magcut = where(zcat.ugriz[2] lt 23.2,ngal)
       hiz = where(zcat[magcut].zbest gt 0.7,nhiz,comp=loz,ncomp=nloz)
       hiz_mostek = get_hiz(zcat.ugriz,/mostek)

; gr vs rz    
       psfile = technotepath+'egs_grz.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.5], yrange=[-0.5,2.5]
       djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
         zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(15), symsize=0.3
       djs_oplot, zcat[hiz].ugriz[2]-zcat[hiz].ugriz[4], $
         zcat[hiz].ugriz[1]-zcat[hiz].ugriz[2], psym=symcat(16), $
         color='orange', symsize=0.3
       im_legend, ['18.5<r<24'], spacing=2.0, /left, /top, box=0, $
         charsize=1.4
       
       rzaxis = range(0.2,1.3,500)
       djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=0, thick=6
       djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=0, thick=6
       djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=0, thick=6
       
       rzcut = 1.2
       blue = where(rzaxis lt rzcut,comp=red)
       djs_oplot, rzaxis[blue], 0.7*rzaxis[blue]+0.05, line=5, thick=6
       djs_oplot, rzaxis[red], 1.4*rzaxis[red]-0.79, line=5, thick=6
       im_legend, ['z<0.7','z>0.7'], /left, /top, box=0, $
         psym=[15,16], color=['','orange']
       
       im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
         line=[0,5], pspacing=1.7, thick=6
       im_plotconfig, psfile=psfile, /psclose, /pdf

stop       
       
; ug vs gr    
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='g - r', ytitle='u - g', xrange=[-0.5,2.0], yrange=[-0.5,3]
       djs_oplot, zcat[loz].ugriz[1,*]-zcat[loz].ugriz[2,*], $
         zcat[loz].ugriz[0,*]-zcat[loz].ugriz[1,*], psym=symcat(15), symsize=0.3
       djs_oplot, zcat[hiz].ugriz[1,*]-zcat[hiz].ugriz[2,*], $
         zcat[hiz].ugriz[0,*]-zcat[hiz].ugriz[1,*], psym=symcat(16), $
         color='orange', symsize=0.3
       im_legend, ['z<0.7','z>0.7'], /left, /top, box=0, $
         psym=[15,16], color=['','orange']
       
; redshift histogram of sources selected using my grz color-cuts 
       zmin = 0.0
       zmax = 1.8
       zbin = 0.1
       nzbin = ceil((zmax-zmin)/zbin)
       zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
       nn = hogg_histogram(zcat.z,[zmin,zmax],nzbin,$
         weight=1.0/zcat.weight)/area
       nn_hiz = hogg_histogram(zcat[colorhiz].z,[zmin,zmax],$
         nzbin,weight=1.0/zcat[colorhiz].weight)/area
       
       im_plotconfig, 6, pos2, height=[6.0,2.5], width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7
       djs_plot, [0], [0], /nodata, position=pos2[*,0], xsty=1, ysty=1, $
         xtitle='', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
         xrange=[zmin,zmax], yrange=[0.0,5000], xtickname=replicate(' ',10)
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn,0], psym=10, thick=6, line=5
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn_hiz,0], psym=10, thick=6, line=0, color='red'

       im_legend, ['Parent Sample (N_{eff}='+strtrim(long(total(nn)*area),2)+')',$
         'grz Sample (N_{eff}='+strtrim(long(total(nn_hiz)*area),2)+')'], $
         line=[5,0], pspacing=1.7, color=['','red'], /right, /top, box=0, $
         charsize=1.4, spacing=2.0, margin=0, thick=6

       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,1], xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='Fraction', $
         xrange=[zmin,zmax], yrange=[0.0,1.1], ytickinterval=0.5
       djs_oplot, !x.crange, [1,1], line=1, thick=6
       djs_oplot, zhist, nn_hiz/nn, psym=10, thick=6
       
; [OII] flux vs redshift
       factor = 1D18
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='F([O II]) (10^{-18} erg / s / cm^{2})', $
         xrange=[0.7,1.6], yrange=[0.1,5000], /ylog, $
         title='grz Sample'
       good = where(oii[colorhiz] gt 0 and oiilimit[colorhiz] eq 0)
       djs_oplot, zcat[colorhiz[good]].z, factor*oii[colorhiz[good]], $
         psym=symcat(16), symsize=0.6, color=cgcolor('dark grey')
       lim = where(oii[colorhiz] gt 0 and oiilimit[colorhiz] eq 1)
       djs_oplot, zcat[colorhiz[lim]].z, factor*oii[colorhiz[lim]], $
         psym=symcat(11,thick=4), symsize=0.5, color=cgcolor('firebrick')
       none = where(oii[colorhiz] lt 0,nnone)
       djs_oplot, zcat[colorhiz[none]].z, replicate(10.0,nnone), $
         psym=symcat(6), symsize=0.8, color=cgcolor('dodger blue')
       djs_oplot, !x.crange, factor*oiicut1*[1,1], line=0, thick=6
       xyouts, 1.45, factor*4.5D-17, textoidl('8\times10^{-17} cgs'), align=0.0, $
         charsize=1.2, /data
       
       im_legend, ['Well-measured','Upper limit (1\sigma)','No measurement'], $
         /right, /bottom, box=0, psym=[16,11,6], $
         color=['dark grey','firebrick','dodger blue'], $
         charsize=1.6
       
; number counts
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
         xrange=[brightcut,faintcut], yrange=alog10([10,5E4]) ;, /ylog
       djs_oplot, magaxis, alog10(total(dndm,/cumu)), line=0, thick=6, $
         psym=10                ;, psym=symcat(16)
       
       ww = where(total(dndm_colorhiz,/cumu) gt 0)
       djs_oplot, magaxis[ww], alog10((total(dndm_colorhiz,/cumu))[ww]), $
         color='orange', line=3, thick=6, psym=10 ;, psym=symcat(15)
       ww = where(total(dndm_oiibright,/cumu) gt 0)
       djs_oplot, magaxis[ww], alog10((total(dndm_oiibright,/cumu))[ww]), color='red', $
         line=1, thick=6, psym=10 ; psym=symcat(6,thick=4), symsize=1.3
       ww = where(total(dndm_oiibright_cor,/cumu) gt 0)
       djs_oplot, magaxis[ww], alog10((total(dndm_oiibright_cor,/cumu))[ww]), color='red', $
         line=5, thick=8, psym=10 ; psym=symcat(6,thick=4), symsize=1.3
       
       rcut = 23.0
;      numcut = 2300.0
;      rcut = interpol(magaxis,total(dndm_oiibright_cor,/cumu),numcut)
       numcut = long(interpol(total(dndm_oiibright_cor,/cumu),magaxis,rcut))
       djs_oplot, [!x.crange[0],rcut], alog10(numcut)*[1,1], line=1, thick=6
       djs_oplot, rcut*[1,1], [!y.crange[0],alog10(numcut)], line=1, thick=6

       xyouts, 18.3, alog10(numcut*1.2), textoidl(string(numcut,$
         format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.4
       
       im_legend, ['Parent Sample','grz Sample','grz && F([OII])>8\times10^{-17} cgs'], $
         color=['','orange','red'], /left, /top, box=0, thick=6, $
         charsize=1.4, spacing=2.0, margin=0, line=[0,3,5], pspacing=1.5
       
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

; ###########################################################################
; build some QAplots of the parent samples for the technote 
    if keyword_set(qaplots_parent) then begin
       allphot = mrdfits(deep2_path(/cat)+'deep2.pcat_ext.fits.gz',1)    
       egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
       allphot = allphot[egs]

       allzcat = read_deep2_zcat(/all)
       egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',negs)
       allzcat = allzcat[egs]
       
       phot = mrdfits(outpath+'egs_parent_phot.fits.gz',1)
       zcat = mrdfits(outpath+'egs_parent_zcat.fits.gz',1)
       zcat_q34 = mrdfits(outpath+'egs_parent_zcat.Q34.fits.gz',1)

       splog, 'Temporary HACK!! Use the maskweight!'
       weight = zcat.mask_weight
       weight_q34 = zcat_q34.mask_weight
;      weight = zcat.targ_weight
;      weight_q34 = zcat_q34.targ_weight

       nphot = n_elements(phot)
       nzcat = n_elements(zcat)
       nzcat_weighted = long(total(1.0/weight))

       nzcat_q34 = n_elements(zcat_q34)
       nzcat_q34_weighted = long(total(1.0/weight_q34))
       splog, nzcat, nzcat_weighted, nzcat_q34, nzcat_q34_weighted

; ra, dec       
       psfile = technotepath+'egs_radec.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7;, ypage=6.5, xpage=8.5

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='R. A. (J2000, deg)', ytitle='Dec (J2000, deg)', $
         xrange=[216.5,213.0], yrange=[51.8,53.8], xtickinterval=1.0;, $
;        title='DEEP2/EGS (Area='+string(area,format='(F5.3)')+' deg^{2})'
       djs_oplot, allphot.ra_deep, allphot.dec_deep, psym=3, color='grey'
       djs_oplot, phot.ra, phot.dec, psym=6, symsize=0.15, color='red'
       djs_oplot, allzcat.ra, allzcat.dec, psym=7, symsize=0.1, color='cyan'
       djs_oplot, zcat.ra, zcat.dec, psym=symcat(16), symsize=0.1, color='blue'

       im_legend, ['PhotParent (N='+strtrim(nphot,2)+')',$
         'zcatParent (N='+strtrim(nzcat_q34,2)+')'],$
         psym=[6,16], color=['red','blue'], /left, /bottom, box=0, $
         charsize=1.4, spacing=2.0, margin=0

       im_legend, ['Matthews+13 (all photometry)','Newman+13 (all spectroscopy)'], $
         psym=[6,6], color=['grey','cyan'], symsize=0.5, /right, /top, box=0, $
         charsize=1.2, spacing=2.0, margin=0

;       im_legend, ['Matthews+13 (N='+strtrim(nphot,2)+')'], $
;         /left, /top, box=0, psym=6, color='red', $
;         position=[0.58,0.91], /norm, charsize=1.4

;      im_legend, ['18.5<r<24','w>0','P_{gal}>0.5'], spacing=2.0, $
;        im_legend, ['18.5<r<24','P_{gal}>0.5'], spacing=2.0, $
;          /left, /top, box=0, position=[0.62,0.87], /norm, charsize=1.4

;        im_legend, ['DEEP2 (N='+strtrim(nspec,2)+')'], $
;          /left, /bottom, box=0, psym=16, color='blue', $
;          position=[0.18,0.35], /norm, charsize=1.4
;        im_legend, ['18.5<r<24'], spacing=1.9, $
;          /left, /bottom, box=0, position=[0.22,0.28], /norm, charsize=1.4
;      im_legend, ['18.5<r<24','Q>=3'], spacing=1.9, $
;        /left, /bottom, box=0, position=[0.22,0.26], /norm, charsize=1.4
;      im_legend, ['18.5<r<24','w>0','Q>=3'], spacing=1.9, $
;        /left, /bottom, box=0, position=[0.22,0.22], /norm, charsize=1.4
       im_plotconfig, psfile=psfile, /psclose, /pdf

; differential number counts 
       dndm = get_dndm(phot.ugriz[2,*],magaxis=magaxis,faintcut=faintcut,$
         brightcut=brightcut)
       zcat_dndm_noweight = get_dndm(zcat_q34.ugriz[2,*],$
         faintcut=faintcut,brightcut=brightcut)
       zcat_dndm = get_dndm(zcat_q34.ugriz[2,*],weight=1.0/weight_q34,$
         faintcut=faintcut,brightcut=brightcut)
       
       psfile = technotepath+'egs_dndm.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
         xrange=[brightcut,faintcut], yrange=alog10([10,1E5]);, /ylog
       djs_oplot, magaxis, alog10(total(dndm,/cumu)), psym=symcat(16)
       djs_oplot, magaxis, alog10(total(zcat_dndm_noweight,/cumu)), $
         color='orange', psym=symcat(15)
       djs_oplot, magaxis, alog10(total(zcat_dndm,/cumu)), color='red', $
         psym=symcat(6,thick=4), symsize=1.3

       im_legend, ['PhotParent (N='+strtrim(nphot,2)+')',$
         'zcatParent (observed, N='+strtrim(nzcat_q34,2)+')',$
         'zcatParent (weighted, N_{eff}='+strtrim(nzcat_q34_weighted,2)+')'], $
         psym=[16,15,6], color=['','orange','red'], /left, /top, box=0, $
         charsize=1.4, spacing=2.0, margin=0
;      im_legend, ['Matthews+13 (N='+strtrim(nphot,2)+')',$
;        'DEEP2 (unweighted, N='+strtrim(nspec,2)+')',$
;        'DEEP2 (weighted, N_{eff}='+strtrim(nspec_weighted,2)+')'], $
;        psym=[16,15,6], color=['','orange','red'], /left, /top, box=0, $
;        charsize=1.4, spacing=2.0, margin=0
       im_plotconfig, psfile=psfile, /psclose, /pdf

; redshift histogram
       psfile = technotepath+'egs_zhist.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       zmin = 0.0
       zmax = 1.8
       zbin = 0.1
       nzbin = ceil((zmax-zmin)/zbin)
       zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
       nn = hogg_histogram(zcat_q34.zbest,[zmin,zmax],nzbin)/area
       nn_cor = hogg_histogram(zcat_q34.zbest,[zmin,zmax],nzbin,$
         weight=1.0/weight_q34)/area
  
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
         xrange=[0.0,2.0], yrange=[0.0,5000]
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn_cor,0], psym=10, thick=6, line=0, color='red'
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn,0], psym=10, thick=6, line=5
;      xyouts, 0.7, 2400.0, 'x15', charsize=1.5, align=0.5, /data
       
       im_legend, ['zcatParent (unweighted, N='+strtrim(nzcat_q34,2)+')',$
         'zcatParent (weighted, N_{eff}='+strtrim(nzcat_q34_weighted,2)+')'], $
         line=[5,0], pspacing=1.7, color=['','red'], /right, /top, box=0, $
         charsize=1.4, spacing=2.0, margin=0       
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif
    
return
end
    
