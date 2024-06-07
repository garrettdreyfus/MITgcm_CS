from graph import folderMap, folderMapRefresh, timeSeriesDashboard, folderMapTimeSeries
import graph
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
                "varysf16":{"specialstring":['sf1','sf5','sf10','sf50','sf100','sf150','highres-sf10'], "marker":["v","v","v","v"] ,"color":["pink","red","orange","black","purple","gray","cyan"],"description":["Different shelf depths"]},\
                #"varysf16":{"specialstring":['highres-sf10','sf10'], "marker":["v","v","v","v"] ,"color":["pink","red","orange","black"],"description":["Different shelf depths"]},\
                #"colderboring16":{"specialstring":['d601'], "marker":["o","o","o"] ,"color":["blue"],"description":["Different shelf depths"]},\
                #"coldernewdens16":{"specialstring":['d601'], "marker":["|","|","|"] ,"color":["green"],"description":["Different shelf depths"]}
                }
#fastExplosionCheck(runsdict)
#folderMapRefresh(runsdict)
folderMapTimeSeries(runsdict)
#folderMap(runsdict)
#fig, axises = plt.subplots(2,3)
#timeSeriesDashboard(runsdict,"",fig,axises)
#folderMapGeneric(gprimeWidth,runsdict)
plt.show()
#plt.show()

#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1d300/results/","Reference",quant="THETA",dim="meridional")
#graph.TSAnim("/jbod/gdf/MITgcm_CS/experiments/coldernewdens16/d601/results/","Reference")
#graph.TSheatmap("/jbod/gdf/MITgcm_CS/experiments/coldernewdens16/d601/results/","Reference")
#graph.TSheatmap("/jbod/gdf/MITgcm_CS/experiments/colderboring16/d601/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/colderboring16/d500/results/","Reference")
#graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/coldernewdens16/d500/results/","Reference")
#plt.show()
#graph.circulationFigure("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
#plt.show()

#folderMapGeneric(steadyStateAverageSimple,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

