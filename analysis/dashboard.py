from graph import folderMap, folderMapRefresh, timeSeriesDashboard, folderMapTimeSeries
import graph
import analysis
import matplotlib.pyplot as plt


#generateRunsTable(runsdict)
#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")
#runsdict = {\
                #"varysf16":{"specialstring":['sf50','sf100','sf150'], "marker":["v","v","v"] ,"color":["red","blue","orange"],"description":["Different shelf depths"]},\
                #"colderboring16":{"specialstring":['d500','d601','d700'], "marker":["o","o","o"] ,"color":["red","blue","orange"],"description":["Different shelf depths"]},\
                #"coldernewdens16":{"specialstring":['d500','d601','d700'], "marker":["|","|","|"] ,"color":["red","blue","orange"],"description":["Different shelf depths"]},\
                #"shelfdepth-16-fix16":{"specialstring":['tcld500','tcld600'], "marker":["o","o","o"] ,"color":["purple","gray","green"],"description":["Different shelf depths"]},\
                #"shelfdepth-16-sf16":{"specialstring":['d500','d600','d700'], "marker":["s","s","s"] ,"color":["gray","black","pink"],"description":["Different shelf depths"]},\
                #"sfwarm16":{"specialstring":['d500','d600','d700'], "marker":["o","o","o"] ,"color":["black","pink","green"],"description":["Different shelf depths"]},\
                #}

runsdict = {\
                #"varysf16":{"specialstring":["sf10","sf3"], "marker":["v","v","v","v"] ,"color":["green","purple","orange","black","purple","gray","cyan"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['imroved-coldstart-warm-sf5','imroved-warmstart-warm-sf5','imroved-warmstart-warm-sf8','imroved-coldstart-warm-sf8','imroved-warmstart-warm-sf10','imroved-coldstart-warm-sf10'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['imroved-coldstart-warm-sf5','imroved-warmstart-warm-sf5','imroved-warmstart-warm-sf8'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['unsalted150','salted150'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['sf10salted250'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                "width16":{"specialstring":['salted75','salted150','salted250','saltedsf25w75','saltedsf25w150','saltedsf25w250','nobathsaltedsf10250','nobathsaltedsf1075','nobathsaltedsf10150'], "marker":["x"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['salted75','salted250','salted150','sf10salted250','nobathsaltedsf10150','nobathsaltedsf1075','nobathsaltedsf10250','saltedsf25w250','saltedsf25w75'], "marker":["x"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","teal","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['saltedsf25w75','saltedsf25w150','saltedsf25w250'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna","orange","black","light green"]*2,"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['saltedsf25w75','saltedsf25w150','saltedsf25w250'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna","orange","black","light green"]*2,"description":["Different shelf depths"]},\
                #"width16":{"specialstring":['salted75','salted150','salted250',], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"icefront16":{"specialstring":['sf2point5front150','sf2point5front75','coldersf50front100'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                "icefront16":{"specialstring":['coldersf50front75','coldersf50front100','coldersf50front150','coldersf4front75','coldersf4front100','coldersf4front150','sf2point5front150','sf2point5front75','sf2point5front100','sf4front150','sf4front75','sf4front100'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"icefront16":{"specialstring":['coldersf50front75','coldersf50front100','coldersf50front150'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"icefront16":{"specialstring":['coldersf50front100'], "marker":["x"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['sf50'], "marker":["v"]*20 ,"color":["tan","green","gray","olive","cyan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['sf3d900','sf3cd700'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['sf3d900','sf3cd700','sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf150','imroved-coldstart-warm-sf8','imroved-warmstart-warm-sf10','imroved-coldstart-warm-sf10'], "marker":["$sf$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                #'imroved-warmstart-warm-sf8',
                "varysf16":{"specialstring":['sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150','sf3d300','imroved-coldstart-warm-sf8','imroved-warmstart-warm-sf8','imroved-warmstart-warm-sf5','imroved-coldstart-warm-sf5','sf3d900','imroved-warmstart-warm-sf10','imroved-coldstart-warm-sf10','sf3cd150','sf3cd450','sf3cd700'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['sf1','sf2','imroved-coldstart-warm-sf5','imroved-warmstart-warm-sf5','imroved-warmstart-warm-sf8'], "marker":["v"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['coldstart-warm-sf20','warm-sf10','sf10'], "marker":["v","v","v","v"] ,"color":["pink","red","orange","black"],"description":["Different shelf depths"]},\
                #"colderboring16":{"specialstring":['d601'], "marker":["o","o","o"] ,"color":["blue"],"description":["Different shelf depths"]},\
                #"coldernewdens16":{"specialstring":['d601'], "marker":["|","|","|"] ,"color":["green"],"description":["Different shelf depths"]}
                }
#fastExplosionCheck(runsdict)

#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3d300/results","Reference",quant="SALT",dim="meridional",selval=150*1000)
#plt.show()

#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/width16/unsalted250/results","reference",200*10**3,quant="SALT",dim="meridional",show=True)
#folderMapRefresh(runsdict)
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results/","sf1")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf2/results/","sf2")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/icefront16/coldersf50front150/results/","coldersf50front150")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/imroved-warmstart-warm-sf8/results/","sf4")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results/","sf5")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3d300/results/","sf3d300")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf50/results/","sf50")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf150/results/","sf150")
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf100/results/","sf100")
#folderMap(runsdict)
#plt.show()
analysis.saltBoxes("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results")
#folderMapTimeSeries(runsdict,"")
#graph.folderMapGeneric(graph.steadyStateAverageSimple,runsdict)
#plt.show()
#plt.show()
#plt.show()
#
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results","reference",200*10**3,quant="RHOAnoma",dim="meridional",show=True,fixcb=False)
#plt.show()
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf50/results/","Reference",quant="THETA",dim="meridional")
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results/","Reference",quant="DENS",dim="meridional")
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results/","Reference",quant="SALT",dim="meridional")
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results/","Reference",quant="SALT",dim="meridional")
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/icefront16/sf4front75/results/","Reference",quant="THETA",dim="meridional")
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results/","Reference",quant="THETA",dim="meridional")
#
#graph.bottomAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results/","Reference",quant="THETA")
#$graph.bottomAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results/","Reference",quant="THETA")
#graph.bottomAnim("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results/","Reference",quant="THETA")

#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/icefront16/coldersf4front100/results","Reference",quant="THETA",dim="meridional")
#plt.show()
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/width16/salted75/results","Reference",quant="THETA",dim="meridional")
#plt.show()
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results","Reference",quant="SALT",dim="meridional")
#plt.show()
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf50/results","Reference",quant="SALT",dim="meridional")
#plt.show()
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results","Reference",quant="THETA",dim="zonal")
#plt.show()

#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3cd150/results","Reference",quant="SALT",dim="meridional")
#plt.show()
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3cd700/results","Reference",quant="SALT",dim="meridional")
#plt.show()
#graph.TSAnim("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results/","Reference")
#graph.TSAnim("/jbod/gdf/MITgcm_CS/experiments/width16/salted150/results/","Reference")
#graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3cd700/results/","High slope")
#graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3cd150/results/","Low slope")
#graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results/","Low salt flux")
#graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results/","High salt flux")
##graph.TSheatmap("/jbod/gdf/MITgcm_CS/experiments/colderboring16/d601/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/nobathsaltedsf10150/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/nobathsalted250/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results/","Reference")

##graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/salted75/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/sf10salted250/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results/","Reference")
#plt.show()
#graph.circulationFigure("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results","Reference")
#graph.circulationFigure("/jbod/gdf/MITgcm_CS/experiments/width16/salted150/results","Reference")
#plt.show()

#folderMapGeneric(steadyStateAverageSimple,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

