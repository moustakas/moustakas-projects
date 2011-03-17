pro compare_kewley_grids
; jm03aug26uofa
; test the Kewley et al. (2001) grids using HII region observations 

    plotsym, 0, 0.2, /fill
    
    regions = read_hii_regions()
    grids = read_kewley_grids()

; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=1
    
    good = where((regions.nii_h_alpha gt -900.0) and (regions.oiii_h_beta gt -900.0))
    xregion = regions[good].nii_h_alpha
    yregion = regions[good].oiii_h_beta

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=2
    
    good = where((regions.sii_h_alpha gt -900.0) and (regions.oiii_h_beta gt -900.0))
    xregion = regions[good].sii_h_alpha
    yregion = regions[good].oiii_h_beta

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=3
    
;   good = where((regions.oi_h_alpha gt -900.0) and (regions.oiii_h_beta gt -900.0))
;   xregion = regions[good].oi_h_alpha
;   yregion = regions[good].oiii_h_beta
;
;   djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=4
    
;   good = where((regions.oii_h_alpha gt -900.0) and (regions.oiii_oii gt -900.0))
;   xregion = regions[good].oiii_oii
;   yregion = regions[good].oii_h_alpha
;
;   djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=5
    
;   good = where((regions.oii_h_alpha gt -900.0) and (regions.nii_sii gt -900.0))
;   xregion = regions[good].nii_sii
;   yregion = regions[good].oii_h_alpha
;
;   djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=6
    
    good = where((regions.nii_6584_oii gt -900.0) and (regions.oiii_5007_oii gt -900.0))
    xregion = regions[good].nii_6584_oii
    yregion = regions[good].oiii_5007_oii

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=7
    
    good = where((regions.nii_oii gt -900.0) and (regions.oiii_h_beta gt -900.0))
    xregion = regions[good].nii_oii
    yregion = regions[good].oiii_h_beta

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=8
    
    good = where((regions.nii_oii gt -900.0) and (regions.R23 gt -900.0))
    xregion = regions[good].R23
    yregion = regions[good].nii_oii

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=9
    
    good = where((regions.oii_h_beta gt -900.0) and (regions.oiii_h_beta gt -900.0))
    xregion = regions[good].oii_h_beta
    yregion = regions[good].oiii_h_beta

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    
    plot_kewley_grids, plot=10
    
    good = where((regions.R23 gt -900.0) and (regions.oiii_oii gt -900.0))
    xregion = regions[good].R23
    yregion = regions[good].oiii_oii

    djs_oplot, xregion, yregion, ps=8
    cc = get_kbrd(1)
; ---------------------------------------------------------------------------    

return
end
