from graph import folderMap, folderMapRefresh, timeSeriesDashboard, folderMapTimeSeries, folderMapCombined
import graph
import analysis
import utils
import matplotlib.pyplot as plt


#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")
runsdict = {\
                  "morediagsenhancedgprime16":{"specialstring":['s-20sf10','s20sf10','s40sf10'], "marker":["$g'$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                  "width16":{"specialstring":['w75sf5','w150sf5','w250sf5','w75sf25','w150sf25','w250sf25'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
                  "varysf16":{"specialstring":['sf3sd300','sf3sd900','sf3cd150','sf3cd450','sf3cd700','sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150'],\
                              "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*9+['$warm$']*7 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*3,"description":["Different shelf depths"]},\
                 }

# officialrunsdict = {\
#                   "morediagsenhancedgprime16":{"specialstring":['s-20sf10','s20sf10','s40sf10'], "marker":["$g'$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
#                   "width16":{"specialstring":['w75sf5','w150sf5','w250sf5','w75sf25','w150sf25','w250sf25','nobathsaltedsf10150','nobathsaltedsf10150','nobathsaltedsf10250'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
#                 #"icefront16":{"specialstring":['coldersf50front75','coldersf50front100','coldersf50front150'], "marker":["$if$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
#                   "varysf16":{"specialstring":['sf3sd300','sf3sd900','sf3cd150','sf3cd450','sf3cd700','sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150','sf1sd300','sf1sd900','warm-sf10','warm-sf3'], "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*3,"description":["Different shelf depths"]},\
#                  }
warmrunsdict = {\
                  "varysf16":{"specialstring":["imroved-coldstart-warm-sf0","imroved-coldstart-warm-sf10","imroved-coldstart-warm-sf5","imroved-coldstart-warm-sf8","imroved-warmstart-warm-sf10","imroved-warmstart-warm-sf5","imroved-warmstart-warm-sf8","sf5","sf8","sf10"],\
                "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*20 ,\
                "color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*3,\
                "description":["Different shelf depths"]},\
                 }
# # runsdict['varysf16']["specialstring"] =  runsdict["varysf16"]["specialstring"] +["imroved-coldstart-warm-sf0","imroved-coldstart-warm-sf10","imroved-coldstart-warm-sf5","imroved-coldstart-warm-sf8","imroved-warmstart-warm-sf10","imroved-warmstart-warm-sf5","imroved-warmstart-warm-sf8"]

# runsdict = warmrunsdict
# graph.folderMapGeneric(graph.connectionPlot,runsdict, ylabel = r'Connectedness Fraction', xlabel = r"$B_{\text{total}} (\frac{m^4}{s^3})$")
# plt.grid(True)
# plt.savefig('connectedness.svg')
# exit()
# runsdict = {\
    # "varysf16":{"specialstring":['sf1','sf2','sf3','sf4','sf5','sf10','sf50','sf100','sf150'], "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                 # }
# runsdict = {\
#                   "morediagsenhancedgprime16":{"specialstring":[], "marker":["$g'$"]*20 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
                   # "width16":{"specialstring":['nobathsaltedsf10150'], "marker":["$w$"]*20 ,"color":["pink","gray","red","orange","black","purple","cyan","green","olive","tan","rosybrown","sienna"],"description":["Different shelf depths"]},\
#                   "varysf16":{"specialstring":[], "marker":["$sd$","$sd$","$cd$","$cd$","$cd$"]+["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
             # }

sfcomprunsdict = {\
                 "varysf16":{"specialstring":['sf10','sf50','sf100','sf150'], "marker":["$sf$"]*10 ,"color":["pink","red","orange","black","purple","cyan","green","gray","olive","tan","rosybrown","sienna"]*2,"description":["Different shelf depths"]},\
}




# graph.overturning_plot("/data/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results","sf1")

# graph.overturning_plot("/data/jbod/gdf/MITgcm_CS/experiments/varysf16/sf10/results","sf10")
# plt.show()
#
#
# # graph.crossSectionAverage("/data/jbod/gdf/MITgcm_CS/experiments/varysf16/sf100/results","sf100",quant="DENS",dim="zonal",selval=155*1000,fixcb=False,savepath='/data/jbod/gdf/MITgcm_CS/analysis/out.png')
# exit()
#utils.generateRunsTable(runsdict)
#exit()
#fastExplosionCheck(runsdict)

# THESE FIGURES ARE IN PAPER
# analysis.letGgoCrazy("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results",0,"sf1")
# graph.buildPortfolio("/data/jbod/gdf/MITgcm_CS/experiments/varysf16/imroved-coldstart-warm-sf0/results","coldstart-warm-sf0")
#graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/width16/w75sf5/results","s40sf10")

# folderMapRefresh(warmrunsdict)
# graph.folderMapGeneric(graph.connectionPlot,runsdict)
# plt.savefig("out.png")
# exit()
# # folderMapRefresh(runsdict)
# folderMap(runsdict,savepath="out.png")
# plt.show()
# exit()

# graph.folderMapGeneric(analysis.breakdown,runsdict,savepath="/data/jbod/gdf/MITgcm_CS/pics/out.png",threed=False,axx=2,axy=2)
# plt.show()
# folderMap(runsdict)
# plt.savefig("finalfit.svg")
# plt.show()
# folderMapTimeSeries(warmrunsdict,"")
# plt.savefig("dashboard.png")
# exit()

# graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
# graph.buildPortfolio("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
# plt.show()

# graph.overturning_plot("/jbod/gdf/MITgcm_CS/experiments/morediagsenhancedgprime16/s40sf10/results","s40sf10")
# plt.show()

# graph.folderMapGeneric(graph.steadyStateAverageSimple,runsdict)
# plt.savefig("out.png")
# exit()

# analysis.overturning("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf1/results")
# plt.show()

#folderMapTimeSeries(runsdict,"")

# graph.folderMapGeneric(graph.steadyStateAverageSimple,runsdict)
# plt.savefig('out.pdf')
#plt.show()
#plt.show()
#
#graph.crossSectionAverage("/jbod/gdf/MITgcm_CS/experiments/morediagvarysf16/sf10/results","reference",200*10**3,quant="RHOAnoma",dim="meridional",show=True,fixcb=False)
#plt.show()
#graph.crossSectionAnim("/jbod/gdf/MITgcm_CS/experiments/varysf16/sf50/results/","Reference",quant="THETA",dim="meridional")

##graph.bottomVtop("/jbod/gdf/MITgcm_CS/experiments/width16/salted75/results/","Reference")
#plt.show()
#graph.circulationFigure("/jbod/gdf/MITgcm_CS/experiments/width16/salted250/results","Reference")
#plt.show()
graph.folderMapGeneric(graph.gprimeTheory,sfcomprunsdict,\
                       xlabel = r'$g^{\prime}_{\text{dc}}~\mathrm{(m/s^2)}$',\
                       ylabel = r'$g^{\prime}_{\text{diagnosed}}~\mathrm{(m/s^2)}$')
plt.savefig('gprime.svg')

# # # #################
# #### Plot comparing S(z)/D to Supper
# #################
graph.folderMapGeneric(graph.deltaRhoCompare,sfcomprunsdict,\
                       xlabel = r'$\Delta \rho~\mathrm{(kg~m^{-3})}$',\
                       ylabel = r'$max(\sigma_0)-min(\sigma_0)~\mathrm{(kg~m^{-3})}$')
# plt.xlim(34.3,34.9)
# plt.ylim(34.3,34.9)
# plt.gca().plot([34.3,34.9],[34.3,34.9],linestyle='dashed')
plt.savefig('out.svg')


# #################
#### Plot comparing S(z)/D to Supper
#################
graph.folderMapGeneric(graph.saltLayerDiff,sfcomprunsdict,\
                       xlabel = r'$\overline{\frac{S(z)}{D} }~\mathrm{(g/kg)}$',\
                       ylabel = r'$\overline{S_\text{upper}}~\mathrm{(g/kg)}$')
plt.xlim(34.3,34.9)
plt.ylim(34.3,34.9)
plt.gca().plot([34.3,34.9],[34.3,34.9],linestyle='dashed')
plt.savefig('sgade.svg')
exit()
