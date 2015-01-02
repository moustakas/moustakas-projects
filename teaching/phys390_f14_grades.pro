pro phys390_f14_grades, alldata, test=test, sendit=sendit, final=final
; jm14sep11siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/390-F14/grades/'
;   path = getenv('TEACHING_DIR')+'/390-F14/grades/'
    
    date = '14dec15' ; update this
    semester = 'Fall 2014'
    class = 'Astronomy 390 - Principles of Astrophysics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys390_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Quizzes','Talk','Final Problem Set']
    weight = [0.30,0.40,0.15,0.15]
    droplowest = [0,1,0,0]

;   keep = where(strmatch(data.last_name,'*Mej*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
