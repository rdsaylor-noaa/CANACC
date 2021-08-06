#!/usr/bin/env pythonw
#
# pdiag.py - interactive diagnostic plot to show simultaneous canopy profiles
#            with SEB modeled and measured fluxes
#
import os
import sys
from libaccess import pltutils
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import seaborn as sns
from datetime import datetime, timedelta
from matplotlib import rcParams
from matplotlib.widgets import Slider
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

rcParams["mathtext.default"] = "regular"

# main
def main(argv=None):
   if argv is None:
      argv = sys.argv
   
   if len(argv) != 2:
      print("usage: %s SIMNAME" % os.path.basename(sys.argv[0]))
      return 2

   # process arguments
   simname = argv[1]

   tslice=0

   dirname = "seb"
   varnames = ["rnet", "h", "le", "g"]
   varlabels = ["R$_{n}$", "H", "LE", "G"]
   scolors = ["darkorange", "green", "royalblue", "peru"]
   outfn = "peflux"
   varunits = "W m$^{-2}$"
   plttitle = "SEB Fluxes"

   # read elapsed hour/datetime key file
   dts, hrs = pltutils.timekeys(simname)

   # adjust model times to match observations
   mdts = []
   pdts = []
   for dt in dts:
       pdts.append(dt)
       mdt = datetime.strftime(dt, "%Y-%m-%d %H:%M:%S")
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

   sebvals = {}
   nhr = 0
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
            sebvals[nhr] = [orn, ole, oh]
         else:
            sebvals[nhr] = [-999., -999., -999.]
         nhr=nhr+1

   fig = plt.figure(figsize=(12, 12))

   # create the SEB plot
   ax = plt.subplot2grid((2,4), (0,0), colspan=4)

   # plot model results
   for varname in varnames:
      plt.plot(dts, varx[varname], color=clrs[varname], linestyle="-", linewidth=lnwdth, label=vlbs[varname])

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
   ptseb = plt.title(plttitle+" - "+simname+" - "+mdts[tslice], fontsize=tfsize, y=tyloc)

   # add legend
   if (len(varnames) > 1):
       plt.legend(loc=4, fontsize=lfszlg, bbox_to_anchor=(0.99, 0.65))

   # add vertical line for time slice location
   ymin, ymax = plt.ylim()
   vline, = plt.plot([pdts[tslice], pdts[tslice]], [ymin, ymax], color="k", linestyle="-", marker="None", lw=1.0)
   plt.ylim(ymin=ymin, ymax=ymax)

   # add SEB values for tslice
   xtseb=0.05
   srn = "{:.1f}".format(sebvals[tslice][0])
   sle = "{:.1f}".format(sebvals[tslice][1])
   sh  = "{:.1f}".format(sebvals[tslice][2])
   prn = plt.text(xtseb,1.16, "R$_{n}$ = "+srn, fontsize=lfszlg, transform=ax.transAxes)
   ple = plt.text(xtseb,1.09, "LE = "+sle, fontsize=lfszlg, transform=ax.transAxes)
   ph  = plt.text(xtseb,1.02, "H  = "+sh,  fontsize=lfszlg, transform=ax.transAxes)

   # read column profile data
   z, atair = pltutils.get1Dvar(simname, "met", "tk")
   z, aqh    = pltutils.get1Dvar(simname, "met", "qh")
   z, atlsun = pltutils.get1Dvar(simname, "canopy", "tlsun")
   z, atlshd = pltutils.get1Dvar(simname, "canopy", "tlshd")
   z, agssun = pltutils.get1Dvar(simname, "canopy", "gssun")
   z, agsshd = pltutils.get1Dvar(simname, "canopy", "gsshd")

   # plot air temperature
   ax1 = plt.subplot2grid((2,4), (1, 0))
   tair     = atair[:, tslice] - 273.15                       # convert from K to C
   tlsun     = atlsun[:, tslice] - 273.15                      # convert from K to C 
   tlshd     = atlshd[:, tslice] - 273.15                      # convert from K to C 
   ptair, = ax1.plot(tair, z, color="darkblue", linestyle="-", marker="None", lw=lnwdth)
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("$^o$C", fontsize=xfsize, labelpad=tlbpad)
   plt.title("T$_{air}$", fontsize=tfsize)
   tamin = min(tair)
   tamax = max(tair)
   tdif = tlsun-tlshd
   dtx = max(abs(tdif))
   dtx = max(5.0, dtx)
   plt.xlim(xmax=tamax+dtx, xmin=tamin-dtx)

   # plot specific humidity
   ax2 = plt.subplot2grid((2,4), (1, 1))
   qh        = aqh[:, tslice]
   pqh, = ax2.plot(qh, z, color="brown", linestyle="-", marker="None", lw=lnwdth)
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("g/kg", fontsize=xfsize, labelpad=tlbpad)
   plt.title("q$_{h}$", fontsize=tfsize)
   plt.xlim(xmax=24.0, xmin=0.0)
   ax2.axes.yaxis.set_ticklabels([])

   # plot leaf temperature
   ax3 = plt.subplot2grid((2,4), (1, 2))
   ptlsun, = ax3.plot(tlsun, z, color="cyan", linestyle="-", marker="None", lw=lnwdth, label="sun")
   ptlshd, = ax3.plot(tlshd, z, color="blue", linestyle="-", marker="None", lw=lnwdth, label="shd")
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("$^o$C", fontsize=xfsize, labelpad=tlbpad)
   plt.legend(loc=4, fontsize=lfszsm, bbox_to_anchor=(0.50, 0.85))
   plt.title("T$_{leaf}$", fontsize=tfsize)
   plt.xlim(xmax=tamax+dtx, xmin=tamin-dtx)
   ax3.axes.yaxis.set_ticklabels([])
   
   # plot stomatal conductance
   ax4 = plt.subplot2grid((2,4), (1, 3))
   gssun     = agssun[:, tslice]
   gsshd     = agsshd[:, tslice]
   pgssun, = ax4.plot(gssun, z, color="green", linestyle="-", marker="None", lw=lnwdth, label="sun")
   pgsshd, = ax4.plot(gsshd, z, color="olive", linestyle="-", marker="None", lw=lnwdth, label="shd")
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("mol/m$^2$-s", fontsize=xfsize, labelpad=tlbpad)
   plt.legend(loc=4, fontsize=lfszsm, bbox_to_anchor=(0.50, 0.85))
   plt.title("g$_{s}$", fontsize=tfsize)
   plt.xlim(xmax=2.0, xmin=0.0)
   ax4.axes.yaxis.set_ticklabels([])

   # make room for the slider
   plt.subplots_adjust(bottom=0.15)

   # make a slider to control the time slice
   slidercolor = "cyan"
   axtime = plt.axes([0.25, 0.03, 0.55, 0.02], facecolor=slidercolor)
   time_slider = Slider(ax=axtime, label="Time", valmin=0, valmax=nts-1, valstep=1, valinit=0)
 
   # update function for the slider
   def update(val):
       tslice=int(time_slider.val)
       vline.set_xdata([pdts[tslice], pdts[tslice]])
       ptseb.set_text(plttitle+" - "+simname+" - "+mdts[tslice])
       srn = "{:.1f}".format(sebvals[tslice][0])
       sle = "{:.1f}".format(sebvals[tslice][1])
       sh  = "{:.1f}".format(sebvals[tslice][2])
       prn.set_text("R$_{n}$ = "+srn)
       ple.set_text("LE = "+sle)
       ph.set_text("H  = "+sh)
       ptair.set_xdata(atair[:, tslice]-273.15)
       pqh.set_xdata(aqh[:, tslice])
       ptlsun.set_xdata(atlsun[:, tslice]-273.15)
       ptlshd.set_xdata(atlshd[:, tslice]-273.15)
       pgssun.set_xdata(agssun[:, tslice])
       pgsshd.set_xdata(agsshd[:, tslice])
       fig.canvas.draw_idle()

   # register the update function with the slider
   time_slider.on_changed(update)

   # allow slider to be moved by "left" and "right" arrow keys
   def onarrow(event):
       if event.key == "left":
           if (time_slider.val == 0):
               val = 0
           else:
               val = time_slider.val - 1
           time_slider.set_val(val)
       elif event.key == "right":
           if (time_slider.val == (nts-1)):
               val = nts - 1
           else:
               val = time_slider.val + 1
           time_slider.set_val(val)
       else:
           pass

   # bind key event to the arrow key handler
   kid = fig.canvas.mpl_connect("key_press_event", onarrow)

   plt.show()

   return 0

if __name__ == "__main__":
   sys.exit(main())
