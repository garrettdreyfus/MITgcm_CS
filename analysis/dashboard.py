from graph import folderMap, folderMapRefresh, timeSeriesDashboard, folderMapTimeSeries
import graph
import matplotlib.pyplot as plt


#generateRunsTable(runsdict)
#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")
runsdict = {\
                "shelfdepth-16-fix16":{"specialstring":['d500','d600','d700'], "marker":["x","x","x"] ,"color":["red","blue","orange"],"description":["Different shelf depths"]},\
                #"shelfdepth-16-fix16":{"specialstring":['tcld500','tcld600'], "marker":["o","o","o"] ,"color":["purple","gray","green"],"description":["Different shelf depths"]},\
                "shelfdepth-16-sf16":{"specialstring":['d500','d600','d700'], "marker":["s","s","s"] ,"color":["white","black","pink"],"description":["Different shelf depths"]},\
                "sfwarm16":{"specialstring":['d500','d600','d700'], "marker":["o","o","o"] ,"color":["black","pink","green"],"description":["Different shelf depths"]},\
                }

#fastExplosionCheck(runsdict)
#folderMapRefresh(runsdict)
folderMapTimeSeries(runsdict)
#folderMap(runsdict)
#fig, axises = plt.subplots(2,3)
#timeSeriesDashboard(runsdict,"",fig,axises)
#folderMapGeneric(gprimeWidth,runsdict)
#plt.show()
#plt.show()
#crossSectionAverage("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
#plt.show()
#graph.circulationFigure("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
plt.show()

#folderMapGeneric(steadyStateAverageSimple,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

