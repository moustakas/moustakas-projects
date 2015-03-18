pro stellarlocus_decals
; jm15mar13siena - build some stellar locus plots

    common com_stellarlocus_dr1, cat

; get the merged tractor catalog (the only cut on this is
; BRICK_PRIMARY='T') 
    qa_dir = getenv('IM_PROJECTS_DIR')+'/decals/dr1/qaplots/'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/decals/isedfit/dr1/'

; read the photometric catalog
    if (n_elements(cat) eq 0) then begin
       if n_elements(cat) eq 0L then cat = $
         mrdfits(isedfit_dir+'tractor-dr1.fits.gz',1)
    endif

    decals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
    maggies = maggies[0:2,*]
    ivarmaggies = ivarmaggies[0:2,*]
    filterlist = filterlist[0:2]
    mag = maggies2mag(maggies,ivarmaggies=ivarmaggies)

; cut the sample    
    keep = where(strtrim(cat.type,2) eq 'PSF' and cat.decam_saturated eq 'F' and $
      mag[1,*] ge 19.0 and mag[1,*] le 20.5)
    label = ['19<r<20.5']
    maggies = maggies[*,keep]

; compute all the color combinations
    colors = stellarlocus_decals_colors(maggies,filterlist='grz')

; Pickles+98    
    pmaggies = stellarlocus_project_pickles(filterlist,pickles=pickles)
    pcolors = stellarlocus_decals_colors(pmaggies,filterlist='grz')

; Kurucz models
    kmaggies = stellarlocus_project_kurucz(filterlist)
    kcolors = stellarlocus_decals_colors(kmaggies,filterlist='grz')

; Padova isochrones
;   pad2 = stellarlocus_read_padova(2)
;   pad3 = stellarlocus_read_padova(3)
    
; make the plots        
    grrange = [-0.7,2.2]
    rzrange = [-0.7,2.4]

    psfile = qa_dir+'qa_stellarlocus_dr1.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, xmargin=[1.3,0.4], height=6.0, width=6.8
; g-r vs r-z
    stellarlocus_decals_scatterplot, colors.rz, colors.gr, pcolors.rz, pcolors.gr, $
      kcolors.rz, kcolors.gr, xtitle='r-z', ytitle='g-r', xrange=rzrange, $
      yrange=grrange, type=pickles.type, feh=pickles.feh, label=label, $
      position=pos
;   djs_oplot, pad2.g-pad2.r, pad2.u-pad2.g, line=0, thick=5.0, color='orange'
;   djs_oplot, pad3.g-pad3.r, pad3.u-pad3.g, line=5, thick=5.0, color='orange'
    
    im_plotconfig, psfile=psfile, /psclose, /png

stop    
    
return
end
