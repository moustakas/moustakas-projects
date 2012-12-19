pro phys130_f12_grades, alldata, test=test, sendit=sendit, $
  midtermplots=midtermplots
; jm12oct05siena - parse the grades for this class 

    path = '~/Dropbox/teaching/12fall/phys130/grades/'
    
    date = '12dec16' ; update this
    semester = 'Fall 2012'
    class = 'Physics 130 - General Physics I'

    if keyword_set(midtermplots) then begin
       m1 = read_gradefile(path+'phys130_midterm1.csv')
       m2 = read_gradefile(path+'phys130_midterm2.csv')
       m3 = read_gradefile(path+'phys130_midterm3.csv')

       psfile = path+'phys130_midterms.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5, charsize=2.5
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xrange=[20,100], yrange=[0,5], xtitle='Percent Grade', $
         ytitle='Number of Students'

       bin = 1.0
       im_plothist, 100*m1.total/m1.total_points, bin=bin, $
         /overplot, color=im_color('black'), /fill, $
         fcolor=im_color('tan'), thick=6, line=0

       im_plothist, 100*m2.total/m2.total_points, bin=bin, $
         /overplot, color=im_color('navy'), /fill, $
         fcolor=im_color('dodger blue'), thick=6, line=0
       im_plothist, 100*m2.total/m2.total_points, bin=bin, $
         /overplot, color=im_color('navy'), thick=6, line=0

       im_plothist, 100*m3.total/m3.total_points, bin=bin, $
         /overplot, color=im_color('firebrick'), /fill, $
         fcolor=im_color('orange'), thick=6, line=0
       im_plothist, 100*m3.total/m3.total_points, bin=bin, $
         /overplot, color=im_color('firebrick'), thick=6, line=0

;      im_plothist, 100*m1.total/m1.total_points, bin=bin, $
;        /overplot, color=im_color('black'), thick=6, line=1

       im_legend, ['Midterm 1','Midterm 2','Midterm 3'], /left, /top, box=0, $
         color=['tan','dodger blue','orange'], line=0, pspacing=1.9, $
         thick=10, charsize=2
       
       im_plotconfig, psfile=psfile, /psclose, /pdf
       return
    endif
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys130_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Lab','Midterm 1','Midterm 2','Midterm 3','Final Exam']
    weight = [0.2,0.15,0.15,0.15,0.15,0.2]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata
    struct_print, alldata

return
end
