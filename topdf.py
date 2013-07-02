import os
pages = 7
for i in range(pages):
    os.system("psselect -p%d plot.ps plot%d.ps"%(i+1,i+1))
    os.system("epstopdf --autorotate=All plot%s.ps"%(i+1))
os.system('pdftk %s cat output plot.pdf'%(' '.join(['plot%d.pdf'%(i+1) for i in range(pages)])))
os.system('rm %s'%(' '.join(['plot%d.*'%(i+1) for i in range(pages)])))
