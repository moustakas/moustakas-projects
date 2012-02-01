pro talk_12feb_siena, keynote=keynote, noevol=noevol
; jm12jan28ucsd - miscellaneous plots for my 2012 Feb talk at Siena

    mfpath = mf_path()
    mfspath = mf_path(/mfs)
    isedpath = mf_path(/isedfit)
    mzpath = mz_path()
    
    datapath = getenv('IM_RESEARCH_DIR')+'/talks/2012/12feb_siena/'
    if keyword_set(keynote) then talkpath = datapath+'keynote/' else $
      talkpath = datapath

; --------------------------------------------------
; build a Montage of RC3 galaxies (must be run on offshore!)
;   hoggdir = '$IM_ARCHIVE_DIR/projects/hogg_rc3/labeled/'
    hoggdir = '$IM_ARCHIVE_DIR/projects/hogg_rc3/'
    pushd, hoggdir
    rc3 = file_search('*.jpg',count=nall)

    size = '150'
    nrow = '10'
    ngal = 100
;   these = shuffle_indx(nall,num_sub=ngal)

; 155, 381, 322 are bad
    these = [139,43,177,492,74,416,134,477,120,423,199,196,159,367,361,490,180,377,$
      443,405,229,378,460,474,233,480,322,88,499,0,105,438,205,189,440,505,235,156,$
      395,295,308,430,5,437,255,227,36,374,94,174,380,125,297,239,168,274,285,52,$
      394,132,369,343,75,415,269,493,65,366,266,410,126,265,325,136,311,301,242,89,$
      109,157,464,296,219,345,277,382,22,113,178,419,241,8,447,432,212,389,463,263,226,472]
    
    outfile = datapath+'rc3_montage.jpg'
    cmd = 'montage -bordercolor black -borderwidth 1 '+ $
      '-tile '+strtrim(nrow,2)+'x'+strtrim(nrow,2)+' -geometry +0+0 '+$
      '-quality 100 -resize '+size+'x'+size+' '+strjoin(rc3[these],' ')+$
      ' '+outfile
;   splog, cmd
    spawn, cmd, /sh

    popd

stop    
    
return
end
