pro phys140_lab_s13_grades, alldata, test=test, sendit=sendit, final=final
; jm12oct05siena - parse the grades for this class 

    path = getenv('DROPBOX_DIR')+'/Teaching/140-S13/grades/'
    
    date = '13mar13' ; update this
    semester = 'Spring 2013'
    class = 'Physics 140 - General Physics II - Lab'
    
; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys140_lab_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = 'Lab'
    weight = 1.0

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, /lab, $
      test=test, sendit=sendit, alldata=alldata, final=final
    struct_print, alldata

return
end
