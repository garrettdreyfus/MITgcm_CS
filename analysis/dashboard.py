from graph import folderMap, folderMapRefresh, timeSeriesDashboard, folderMapTimeSeries
import graph
import analysis
import utils
import matplotlib.pyplot as plt


#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")
runsdict = {\
                  "morediagsenhancedgprime16":{"specialstring":['s-20sf10','s20sf10','s40sf10'], "marker":["$g'$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                  "width16":{"specialstring":['w75sf5','w150sf5','w250sf5','w75sf25','w150sf25','w250sf25'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                #"icefront16":{"specialstring":['coldersf50front75','coldersf50front100','coldersf50front150'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                  "varysf16":{"specialstring":['sf3sd300','sf3sd900','sf3cd150','sf3cd450','sf3cd700','sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150'], "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                 }
# runsdict = {\
#                   "morediagsenhancedgprime16":{"specialstring":['s-20sf10','s20sf10','s40sf10'], "marker":["$g'$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
#                   "width16":{"specialstring":['w75sf5','w150sf5','w250sf5','w75sf10','w150sf10','w250sf10','w75sf25','w150sf25','w250sf25'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
#                   "icefront16":{"specialstring":['coldersf50front75','coldersf50front100','coldersf50front150'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
#                   "varysf16":{"specialstring":['sf3sd300','sf3sd900','sf3cd150','sf3cd450','sf3cd700','sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150'], "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
#                  }
# 
#runsdict = {\
                 # "varysf16":{"specialstring":['sf5','sf10','sf50','sf100','sf150'], "marker":["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
# }

#utils.generateRunsTable(runsdict)
#exit()
#fastExplosionCheck(runsdict)

# THESE FIGURES ARE IN PAPER
# analysis.letGgoCrazy("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results",0,"sf1")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results","sf1")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results","sf10")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/width16/w250sf5/results","w250sf5")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results","sf10")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/width16/w250sf5/results","sf10")

# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf150/results","sf150")
# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
#graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/width16/w75sf5/results","s40sf10")

# folderMapRefresh(runsdict)
#graph.folderMapGeneric(graph.connectionPlot,runsdict,ylabel="Connectedness Fraction",xlabel=r'$B_{\text{total}}$',savepath="/jbod/gdf/MITgcm_CS/pics/connectionvsB0.svg")
# plt.show()
#folderMapTimeSeries(runsdict,"")

# graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/sf1bump20points/results","Reference",quant="SALT",dim="meridional",selval=150*1000,fixcb=True)
# graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","Reference",show=True)
# graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
# plt.show()

# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
# plt.show()

# graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results","Reference",show=True)
# plt.show()

# graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf100/results","Reference",show=True)
# plt.show()
# graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results","Reference",show=True)
# plt.show()
##graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/sf1bump-20points/results","Reference",quant="SALT",dim="meridional",selval=150*1000,fixcb=True)
#plt.show()
#
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/width16/unsalted250/results","reference",200*10**3,quant="SALT",dim="meridional",show=True)
#graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/sf10bump-20points/results/","sf10bump-20points")
# folderMap(runsdict,savepath="/jbod/gdf/MITgcm_CS/pics/meltevaluation.svg")
# plt.show()
#analysis.saltBoxes("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results")
#analysis.overturning("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results")
#plt.show()
# plt.show()


# analysis.overturning("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results")
# plt.show()

#folderMapTimeSeries(runsdict,"")

graph.folderMapGeneric(graph.steadyStateAverageSimple,runsdict)
plt.show()
#plt.show()
#plt.show()
#
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results","reference",200*10**3,quant="RHOAnoma",dim="meridional",show=True,fixcb=False)
#plt.show()
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf50/results/","Reference",quant="THETA",dim="meridional")

#graph.TSAnim("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results/","Reference")
#graph.volumetricTS("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf3cd700/results/","High slope")
##graph.TSheatmap("/jbod/gdf/MITgcm_CS/experiments/colderboring16/d601/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf5/results/","Reference")

##graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/salted75/results/","Reference")
#plt.show()
#graph.circulationFigure("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results","Reference")
#plt.show()

# graph.folderMapGeneric(graph.connectionPlot,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#plt.show()
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

