#!/usr/bin/env pythonw
#
# peflux.py - plot modeled energy fluxes vs data from Coweeta
#
import os
import sys
from libaccess import pltutils
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import seaborn as sns
from datetime import datetime, timedelta
import csv

# set figure formatting parameters
tfsize   = 18     # plot title font size
tyloc    = 1.02   # title y location
lfszsm   = 10     # legend font size small
lfszlg   = 14     # legend font size large
yfsize   = 18     # y-axis title font size
ylabpad  = 10     # y-axis title padding
xfsize   = 18     # x-axis title font size
xlabpad  = 10     # x-axis title padding
tlmaj    = 6      # major tick length
tlmin    = 4      # minor tick length
tlbsize  = 17     # tick label font size
tlbpad   = 3      # tick label padding
lnwdth   = 2.0    # linewidth
msize    = 4      # marker size

# main
def main(argv=None):
   if argv is None:
      argv = sys.argv
   
   if len(argv) != 3:
      print("usage: %s SIMNAME OUTTYPE" % os.path.basename(sys.argv[0]))
      return 2

   # process arguments
   simname = argv[1]
   outtype = argv[2]

   dirname = "seb"
   varnames = ["rnet", "h", "le", "g"]
   varlabels = ["R$_{n}$", "H", "LE", "G"]
   scolors = ["darkorange", "green", "royalblue", "peru"]
   outfn = "peflux"
   varunits = "W m$^{-2}$"
   plttitle = "Energy Fluxes"

   # read elapsed hour/datetime key file
   dts, hrs = pltutils.timekeys(simname)

   # adjust model times to match observations
   mdts = []
   for dt in dts:
#      mdt = dt - timedelta(hours=1.0)
       mdt = dt
       mdts.append(mdt)

   # get model values to plot
   varx = {}
   vlbs = {}
   clrs = {}
   for varname, varlabel, scolor in zip(varnames, varlabels, scolors):
       dat = pltutils.get0Dvar(simname, dirname, varname)
       varx[varname] = dat
       clrs[varname] = scolor
       vlbs[varname] = varlabel

   nts = len(dts)
   startdt = dts[0]
   enddt = dts[nts-1]

   print(startdt)
   print(enddt)

   # get observations to plot
   obsfn = "./obs/Coweeta_Meteo_Hourly_V2.1.csv"
   fh = open(obsfn, "r")
   reader = csv.DictReader(fh)

   odts = []
   orns = []
   oles = []
   ohs  = [] 

   for row in reader:
      year = int(row["year"])
      doy  = int(row["doy"])
      hour = int(row["hour"])
      odt  = datetime(year, 1, 1) + timedelta(days=doy-1, hours=hour)
      if ( (odt >= startdt) and (odt <= enddt) ):
         orn     = float(row["Rn"])
         ole     = float(row["LE"])
         oh      = float(row["H"])
         if (orn != -999. and ole != -999. and oh != -999.):
            odts.append(odt)
            orns.append(orn)
            oles.append(ole)
            ohs.append(oh) 

   # create the plot
   fig, ax = plt.subplots(1, 1, figsize=(12, 6))

   # plot model results
   for varname in varnames:
      plt.plot(mdts, varx[varname], color=clrs[varname], linestyle="-", linewidth=lnwdth, label=vlbs[varname])

   # plot observations
   plt.plot(odts, orns, color=scolors[0], linestyle="None", marker="o", markersize=msize)
   plt.plot(odts, ohs,  color=scolors[1], linestyle="None", marker="o", markersize=msize)
   plt.plot(odts, oles, color=scolors[2], linestyle="None", marker="o", markersize=msize)

   # take care of time formatting on x-axis
   days = mdates.DayLocator()
   ax.xaxis.set_major_locator(days)
   ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
   if (nts > 48):
       hours = mdates.HourLocator(byhour=range(24), interval=4)
   else:
       hours = mdates.HourLocator(byhour=range(24), interval=1)
   ax.xaxis.set_minor_locator(hours)

   # set standard formatting
   pltutils.setstdfmts(ax, tlmaj, tlmin, tlbsize, tlbpad)

   # set y-axis label
   plt.ylabel(varunits, fontsize=yfsize, labelpad=ylabpad)

   # add plot title
   plt.title(plttitle+" - "+simname, fontsize=tfsize, y=tyloc)

   # add legend
   if (len(varnames) > 1):
       plt.legend(loc=4, fontsize=lfszlg, bbox_to_anchor=(0.99, 0.65))

   # create output
   pltutils.pltoutput(simname, outfn, outtype)

   return 0

if __name__ == "__main__":
   sys.exit(main())
