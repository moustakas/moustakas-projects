pro stellarlocus_decals
; jm15mar13siena - build some stellar locus plots

    common com_stellarlocus_dr1, cat

; get the merged tractor catalog (the only cut on this is
; BRICK_PRIMARY='T') 
    dr1dir = getenv('DECALS_DIR')+'/'
    qadir = dr1dir+'qaplots/'
    file_mkdir, qadir
;   qa_dir = getenv('IM_PROJECTS_DIR')+'/decals/dr1/qaplots/'
;   isedfit_dir = getenv('IM_PROJECTS_DIR')+'/decals/isedfit/dr1/'

; read the photometric catalog
    if (n_elements(cat) eq 0) then begin
       allbrick = file_basename(file_search(dr1dir+'tractor/2[3-4]?/',/test_dir))
;      allbrick = allbrick[0:2]
       nbrick = n_elements(allbrick)
       for ii = 0L, nbrick-1 do begin
;         delvarx, cat
          catfile = file_search(dr1dir+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
          for ic = 0L, ncat-1 do begin
             print, format='("Brick ",I0,"/",I0," Cat ",I0,"/",I0, A10,$)', $
               ii+1, nbrick, ic+1, ncat, string(13b)
             cat1 = mrdfits(catfile[ic],1,/silent)
             these = where(cat1.brick_primary eq 'T' and total(cat1.decam_nobs[[1,2,4]] gt 0,1) eq 3)
             if these[0] ne -1 then begin
                cat1 = cat1[these]
;               cat1 = struct_trimtags(cat1[these],except='sdss_*')
                if n_elements(cat) eq 0L then cat = cat1 else cat = [temporary(cat),temporary(cat1)]
             endif
          endfor
;         if n_elements(cat) ne 0L then stop          
       endfor 
       ngal = n_elements(cat)
;      cat = mrdfits(isedfit_dir+'tractor-dr1.fits.gz',1)
    endif

    decals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist, /decam_grz
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

    psfile = qadir+'qa_stellarlocus_dr1.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2, xmargin=[1.3,0.4], height=6.0, width=6.8
; g-r vs r-z
    stellarlocus_decals_scatterplot, colors.rz, colors.gr, pcolors.rz, pcolors.gr, $
      kcolors.rz, kcolors.gr, xtitle='r-z', ytitle='g-r', xrange=rzrange, $
      yrange=grrange, type=pickles.type, feh=pickles.feh, label=label, $
      position=pos
;   djs_oplot, pad2.g-pad2.r, pad2.u-pad2.g, line=0, thick=5.0, color='orange'
;   djs_oplot, pad3.g-pad3.r, pad3.u-pad3.g, line=5, thick=5.0, color='orange'
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
