#!/usr/bin/env pythonw
#
# pmet.py - interactive diagnostic plot to show selected met profiles
#               at a particular hour
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

# colors
colors = ["darkgray", "peru", "red", "darkorange", "royalblue", "green", "violet", "magenta", "cyan", "olive"]

# formatting
tfsize = 20    # plot title font size
tyloc  = 0.95  # title y location
lfsize = 12    # legend font size
yfsize = 18    # y-axis title font size
xfsize = 16    # x-axis title font size
tlmaj  = 6     # major tick length
tlmin  = 4     # minor tick length
tlbsize= 17    # tick label size
msize  = 6     # marker size
lnwdth = 2     # linewidth
xtpad = 0      # x title padding
ytpad = 0      # y title padding
tlbpad = 0     # tick label padding

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

   # read elapsed hour/datetime key file
   dts, hrs = pltutils.timekeys(simname)
   nts = len(dts)
   datetimes = []
   for dt in dts:
       datetimes.append(datetime.strftime(dt, "%Y-%m-%d %H:%M:%S"))

   # read data files

   # cair
   z, acair = pltutils.get1Dvar(simname, "met", "cair")
   cair     = acair[:, tslice]

   # h2o
   z, ah2o  = pltutils.get1Dvar(simname, "met", "h2o")
   h2o      = ah2o[:, tslice]

   # kv
   z, akv   = pltutils.get1Dvar(simname, "met", "kv")
   kv       = akv[:, tslice]

   # pmb
   z, apmb  = pltutils.get1Dvar(simname, "met", "pmb")
   pmb      = apmb[:, tslice]

   # ppfd
   z, appfd = pltutils.get1Dvar(simname, "met", "ppfd")
   ppfd     = appfd[:, tslice]

   # qh
   z, aqh   = pltutils.get1Dvar(simname, "met", "qh")
   qh       = aqh[:, tslice]

   # rh
   z, arh   = pltutils.get1Dvar(simname, "met", "rh")
   rh       = arh[:, tslice]

   # tk
   z, atk   = pltutils.get1Dvar(simname, "met", "tk")
   tk       = atk[:, tslice]

   # ubar
   z, aubar = pltutils.get1Dvar(simname, "met", "ubar")
   ubar     = aubar[:, tslice]

   # create the plots
   fig = plt.figure(figsize=(12,12))

   # air density (molec/cm3)
   ax1 = fig.add_subplot(2,4,1)
   pcair, = ax1.plot(cair, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="c$_{air}$")
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.85,0.05))
#  plt.xlim(xmax=1.0, xmin=0.0)
   plt.ylim(ymax=43.0, ymin=0.0)

   # water density (molecules/cm3)
   ax2 = fig.add_subplot(2,4,2)
   ph2o, = ax2.plot(h2o, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="H$_2$O")
   pltutils.setstdfmts(ax2, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.55,0.05))
#  plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=43.0, ymin=0.0)

   # eddy diffusivity (cm2/s)
   ax3 = fig.add_subplot(2,4,3)
   pkv, = ax3.plot(kv, z, color=colors[4], linestyle="-", marker="None", linewidth=lnwdth, label="K$_{v}$")
   pltutils.setstdfmts(ax3, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
   plt.ylim(ymax=43.0, ymin=0.0)

   # air pressure (mb)
   ax4 = fig.add_subplot(2,4,4)
   ppmb, = ax4.plot(pmb, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="p$_{air}$")
   pltutils.setstdfmts(ax4, tlmaj, tlmin, tlbsize, tlbpad)
#  plt.xlabel("$\mu$mol m$^{-2}$ s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
#  plt.xlim(xmax=20., xmin=-1.)
   plt.ylim(ymax=43.0, ymin=0.0)

   # ppfd
   ax5 = fig.add_subplot(2,4,5)
   pppfd, = ax5.plot(ppfd, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="PPFD")
   pltutils.setstdfmts(ax5, tlmaj, tlmin, tlbsize, tlbpad)
#  plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
#  plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=43.0, ymin=0.0)

   # qh
   ax6 = fig.add_subplot(2,4,6)
   pqh, = ax6.plot(qh, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="q$_h$")
   pltutils.setstdfmts(ax6, tlmaj, tlmin, tlbsize, tlbpad)
#  plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
#  plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=43.0, ymin=0.0)

   # rh
   ax7 = fig.add_subplot(2,4,7)
   prh, = ax7.plot(rh, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="RH")
   pltutils.setstdfmts(ax7, tlmaj, tlmin, tlbsize, tlbpad)
#  plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
#  plt.xlim(xmax=500.0, xmin=-20.0)
   plt.ylim(ymax=43.0, ymin=0.0)

   # ubar
   ax8 = fig.add_subplot(2,4,8)
   pubar, = ax8.plot(ubar, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="u$_{mean}$")
#  xmin, xmax = ax8.get_xlim()
#  plt.xticks(np.arange(xmin, xmax+1, 1000.))
   pltutils.setstdfmts(ax8, tlmaj, tlmin, tlbsize, tlbpad)
#  plt.xlabel("s m$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
#  plt.xlim(xmax=2000., xmin=-100.)
   plt.ylim(ymax=43.0, ymin=0.0)

   psupttle = plt.suptitle(simname+" - "+datetimes[tslice], fontsize=tfsize, y=tyloc)

   # make room for the slider
   plt.subplots_adjust(bottom=0.15)

   # make a slider to control the time slice
   slidercolor="cyan"
   axtime = plt.axes([0.25, 0.03, 0.55, 0.02], facecolor=slidercolor)
   time_slider = Slider(ax=axtime, label="Time", valmin=0, valmax=nts-1, valstep=1, valinit=0)

   # update function for the slider
   def update(val):
       tslice = int(time_slider.val)
       psupttle.set_text(simname+" - "+datetimes[tslice])
       xmax = np.amax(acair[:, tslice])
       xmin = np.amin(acair[:, tslice])
       pcair.set_xdata(acair[:, tslice])
       ax1.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

       xmax = np.amax(ah2o[:, tslice])
       xmin = np.amin(ah2o[:, tslice])
       ph2o.set_xdata(ah2o[:, tslice])
       ax2.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

       xmax = np.amax(akv[:, tslice])
       pkv.set_xdata(akv[:, tslice])
       ax3.set_xlim(xmax=xmax*1.1, xmin=-20.0)

       xmax = np.amax(apmb[:, tslice])
       xmin = np.amin(apmb[:, tslice])
       ppmb.set_xdata(apmb[:, tslice])
       ax4.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)
       
       xmax = np.amax(appfd[:, tslice])
       xmin = np.amin(appfd[:, tslice])
       pppfd.set_xdata(appfd[:, tslice])
       ax5.set_xlim(xmax=2500., xmin=-20.0)

       xmax = np.amax(aqh[:, tslice])
       xmin = np.amin(aqh[:, tslice])
       pqh.set_xdata(aqh[:, tslice])
       ax6.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

       xmax = np.amax(arh[:, tslice])
       xmin = np.amin(arh[:, tslice])
       prh.set_xdata(arh[:, tslice])
       ax7.set_xlim(xmax=100., xmin=0.0)

       xmax = np.amax(aubar[:, tslice])
       xmin = np.amin(aubar[:, tslice])
       pubar.set_xdata(aubar[:, tslice])
       ax8.set_xlim(xmax=xmax*1.1, xmin=0.0)

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
