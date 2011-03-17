pro build_sfrm_tables
; jm10feb16ucsd - build the latex tables for the SFRM paper

    emulateapj = 1
    
    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    mf_fit = mrdfits(sfrmpath+'mf_fit_all.fits.gz',1)
    mf_data = mrdfits(sfrmpath+'mf_data_all.fits.gz',1)
    limits = mrdfits(sfrmpath+'sfrm_limits.fits.gz',1,/silent)
    zbins = sfrm_zbins(nzbins)

; ---------------------------------------------------------------------------    
; mass functions for all, quiescent, and active galaxies in six
; redshift bins

; establish the reference mass    
    refgood = where((mf_data[0].phi gt -999.0) and $
      (mf_data[0].mass gt limits.minmass[0]),nrefmass)
    refmass = mf_data[0].mass[refgood]

    scale = 1E3
    
; build the table    
    colhead1 = '\colhead{'+[ ['$\log(\mathcal{M})$']+'} & ',['$\Phi$']+'} \\ ']
    colunits = '\colhead{'+[ ['$(\mathcal{M}_{\sun})$']+'} & ',$
      ['$(10^{-3}\,{\rm Mpc}^{-3}\,{\rm dex}^{-1}$)']+'} ']
    texcenter = ['c','c']

;   tablenotetext = [$
;     '{a}{Fits were performed using a ``robust'' ordinaryleast-squares bisector fit.}',$
;     '{b}{$I$-band magnitudes are on the Vega system, while $g$- and $r$-band magnitudes are AB-relative.}']
    
; write it out    
    texfile = paperpath+'table_mf.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
       printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{Mass Functions\label{table:mf}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

;   for ii = 0, 0 do begin
    for ii = 0, nzbins-1 do begin
       printf, lun, '\cutinhead{$z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2)+'$}'
       for jj = 0, nrefmass-1 do begin
          phi = mf_data[ii].phi[refgood[jj]]
          phierr = mf_data[ii].phierr[refgood[jj]]
          if (phi gt 0.0) then begin
             strphi = '$'+strtrim(string(scale*phi,format='(F12.3)'),2)+$
               '\pm'+strtrim(string(scale*phierr,format='(F12.3)'),2)+'$'
          endif else strphi = '\nodata'
          strmass = '$'+strtrim(string(refmass[jj],format='(F12.2)'),2)+'$'
          printf, lun, strmass+' & '+strphi+' \\'
       endfor          
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    printf, lun, '\clearpage'
    printf, lun, '\end{landscape}'
    free_lun, lun

stop    
    
return
end
    
