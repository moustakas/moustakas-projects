PRO lhcn_lhacor,pspath,hcndust,lhcnerr,lhcnrange,lharange3,hcnnodust,$
  nametag=nametag,sfrrange2,sfrhcnrange,$
  postscript=postscript,encapsulated=encapsulated,$
  postthick,postthick2

; ---------------------------------------------------------------------------    
; L(HCN) vs L(Ha)_cor
; ---------------------------------------------------------------------------    

psname = 'LHCN.LHacor'
im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.3, encapsulated=encapsulated

pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
xmargin=[1.5,0.4], ymargin=[1.1,1.1], xpage=9.5, ypage=9.0, $
position=pos, /normal

; plotting variables
hcnsfsym = 108 & hcnsfcolor = 'red' & hcnsfpsize = 1.5 & hcnsffill = 1
hcnsfsym_ul = 108 & hcnsfcolor_ul = 'red' & hcnsfpsize_ul = 4.0 
hcnsffill_ul = 1

hcnagnsym = 105 & hcnagncolor = 'blue' & hcnagnpsize = 1.8 & hcnagnfill = 1
hcnagnsym_ul = 105 & hcnagncolor_ul = 'blue' & hcnagnpsize_ul = 4.0 
hcnagnfill_ul = 1

hcnagnsfsym = 106 & hcnagnsfcolor = 'green' & hcnagnsfpsize = 1.8 
hcnagnsffill = 1 & hcnagnsfsym_ul = 106 & hcnagnsfcolor_ul = 'green' 
hcnagnsfpsize_ul = 4.0 & hcnagnsffill_ul = 1

hcnallsym = 108 & hcnallcolor = 'black' & hcnallpsize = 1.0 & hcnallfill = 1
hcnsym = 108 & hcncolor = 'red' & hcnpsize = 1.0 & hcnfill = 1
atlassym = 106 & atlascolor = 'grey' & atlaspsize = 1.0 & atlasfill = 1

indx = where((hcndust.gao_lhcn gt 0.0) and (hcnnodust.sfr_h_alpha gt -900.0),nindx)
indx_ul = where((hcndust.gao_lhcn lt 0.0) and (hcnnodust.sfr_h_alpha gt -900.0),nindx_ul)

x = hcndust[indx].gao_lhcn
xerr = x*0.0 + lhcnerr
y = hcnnodust[indx].h_alpha_lum[0];hcnnodust[indx].sfr_h_alpha
yerr = hcnnodust[indx].h_alpha_lum[1]

x_ul = -hcndust[indx_ul].gao_lhcn
xerr_ul = x_ul*0.0 + lhcnerr
y_ul = hcnnodust[indx_ul].h_alpha_lum[0]
yerr_ul = hcnnodust[indx_ul].h_alpha_lum[1]

agn = where(strmatch(hcndust[indx].bpt_nii_mixture_class,$
  'AGN*',/fold),nagn,comp=sf,ncomp=nsf)
sf = where(strmatch(hcndust[indx].bpt_nii_mixture_class,$
  'HII *',/fold),nsf)
agnsf = where(strmatch(hcndust[indx].bpt_nii_mixture_class,$
  'HII/AGN*',/fold) or strmatch(hcndust[indx].bpt_nii_mixture_class,$
  'Unknown',/fold),nagnsf)
agn_ul = where(strmatch(hcndust[indx_ul].bpt_nii_mixture_class,$
  'AGN*',/fold),nagn_ul,comp=sf_ul,ncomp=nsf_ul)
agnsf_ul = where(strmatch(hcndust[indx_ul].bpt_nii_mixture_class,$
  'HII/AGN*',/fold) or strmatch(hcndust[indx_ul].bpt_nii_mixture_class,$
  'Unknown',/fold),nagnsf_ul)
sf_ul = where(strmatch(hcndust[indx_ul].bpt_nii_mixture_class,$
  'HII *',/fold),nsf_ul)

; make the plot!    

xtitle = 'log L(HCN) [L'+sunsymbol()+']'
ytitle = 'log [L(H\alpha)_{cor}] [L'+sunsymbol()+']'
ytitle2 = textoidl('log \psi[L(H\alpha)_{cor}] [M'+sunsymbol()+' yr^{-1}]')

xrange = lhcnrange
xrange2 = sfrhcnrange
;yrange = sfrrange2

yrange = [6.9,9.6];lharange3
;yrange2 = sfrrange2
yrange2 = sfrrange2

extra = {YSTYLE: 8, XSTYLE: 8}


atlas1d_lineplot, x[sf], y[sf], xerr[sf], yerr[sf], plottype=4, $
  postscript=postscript, $
  xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
  charsize=charsize_8, position=pos[*,0], $;xtickname=replicate(' ',10), $
  atlassym=hcnsfsym, atlaspsize=hcnsfpsize, atlascolor=hcnsfcolor, $
  atlasfill=hcnsffill,_extra=extra
AXIS,YAXIS=1,yrange=yrange2,charsize=2,charthick=postthick,$
  ythick=postthick,$
  ytitle=ytitle2,ystyle=1
atlas1d_lineplot, x[agn], y[agn], xerr[agn], yerr[agn], plottype=4, /overplot, $
  atlassym=hcnagnsym, atlaspsize=hcnagnpsize, atlascolor=hcnagncolor, $
  atlasfill=hcnagnfill
atlas1d_lineplot, x[agnsf], y[agnsf], xerr[agnsf], yerr[agnsf], plottype=4, $
  /overplot, atlassym=hcnagnsfsym, atlaspsize=hcnagnsfpsize, $
  atlascolor=hcnagnsfcolor, atlasfill=hcnagnsffill
AXIS,XAXIS=1,xrange=xrange2,charsize=2,charthick=postthick,xthick=postthick,$
  xtitle=textoidl('log \psi[L(HCN)] [M'+sunsymbol()+' yr^{-1}]'),xstyle=1

; upper limits    

if (nagn_ul ne 0L) then atlas1d_lineplot, x_ul[agn_ul], y_ul[agn_ul], $
  xerr_ul[agn_ul], yerr_ul[agn_ul], plottype=1, /overplot, $
  atlassym=114, atlaspsize=hcnagnpsize_ul, atlascolor=hcnagncolor_ul, $
  atlasfill=hcnagnfill_ul, symthick=postthick2, errthick=postthick2
if (nsf_ul ne 0L) then atlas1d_lineplot, x_ul[sf_ul], y_ul[sf_ul], $
  xerr_ul[sf_ul], yerr_ul[sf_ul], plottype=1, /overplot, $
  atlassym=114, atlaspsize=hcnsfpsize_ul, atlascolor=hcnsfcolor_ul, $
  atlasfill=hcnsffill_ul, symthick=postthick2, errthick=postthick2
if (nagnsf_ul ne 0L) then atlas1d_lineplot, x_ul[agnsf_ul], $
  y_ul[agnsf_ul], xerr_ul[agnsf_ul], yerr_ul[agnsf_ul], plottype=1, /overplot, $
  atlassym=114, atlaspsize=hcnagnsfpsize_ul, atlascolor=hcnagnsfcolor_ul, $
  atlaagnsfill=hcnagnsffill_ul, symthick=postthick2, errthick=postthick2

im_openclose, postscript=postscript, /close

END
