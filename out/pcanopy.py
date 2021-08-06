#!/usr/bin/env pythonw
#
# pcanopy.py - interactive diagnostic plot to show selected canopy profiles
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
lfsize = 10    # legend font size
yfsize = 16    # y-axis title font size
xfsize = 14    # x-axis title font size
tlmaj  = 6     # major tick length
tlmin  = 4     # minor tick length
tlbsize= 14    # tick label size
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

   # fsun
   z, afsun = pltutils.get1Dvar(simname, "canopy", "fsun")
   fsun     = afsun[:, tslice]

   # fshd
   z, afshd = pltutils.get1Dvar(simname, "canopy", "fshd")
   fshd     = afshd[:, tslice]

   # ppfdsun
   z, appfdsun = pltutils.get1Dvar(simname, "canopy", "ppfdsun")
   ppfdsun     = appfdsun[:, tslice]

   # ppfdshd
   z, appfdshd = pltutils.get1Dvar(simname, "canopy", "ppfdshd")
   ppfdshd     = appfdshd[:, tslice]

   # nirsun
   z, anirsun = pltutils.get1Dvar(simname, "canopy", "nirsun")
   nirsun     = anirsun[:, tslice]

   # nirshd
   z, anirshd = pltutils.get1Dvar(simname, "canopy", "nirshd")
   nirshd     = anirshd[:, tslice]

   # lwup
   z, alwup = pltutils.get1Dvar(simname, "canopy", "lwup")
   lwup     = alwup[:, tslice]

   # lwdn
   z, alwdn = pltutils.get1Dvar(simname, "canopy", "lwdn")
   lwdn     = alwdn[:, tslice]

   # rtsun
   z, artsun = pltutils.get1Dvar(simname, "canopy", "rtsun")
   rtsun     = artsun[:, tslice]

   # rtshd
   z, artshd = pltutils.get1Dvar(simname, "canopy", "rtshd")
   rtshd     = artshd[:, tslice]

   # rasun
   z, arasun = pltutils.get1Dvar(simname, "canopy", "rabssun")
   rasun     = arasun[:, tslice]

   # rashd
   z, arashd = pltutils.get1Dvar(simname, "canopy", "rabsshd")
   rashd     = arashd[:, tslice]

   # rssun
   z, arssun = pltutils.get1Dvar(simname, "canopy", "rssun")
   rssun     = arssun[:, tslice]

   # rsshd
   z, arsshd = pltutils.get1Dvar(simname, "canopy", "rsshd")
   rsshd     = arsshd[:, tslice]

   # tlsun
   z, atlsun = pltutils.get1Dvar(simname, "canopy", "tlsun")
   tlsun     = atlsun[:, tslice] - 273.15     # convert from K to C

   # tlshd
   z, atlshd = pltutils.get1Dvar(simname, "canopy", "tlshd")
   tlshd     = atlshd[:, tslice] - 273.15     # convert from K to C

   # gssun
   z, agssun = pltutils.get1Dvar(simname, "canopy", "gssun")
   gssun     = agssun[:, tslice]

   # gsshd
   z, agsshd = pltutils.get1Dvar(simname, "canopy", "gsshd")
   gssun     = agssun[:, tslice]

   # rssun
   z, arssun = pltutils.get1Dvar(simname, "canopy", "rssun")
   rssun     = arssun[:, tslice]

   # rsshd
   z, arsshd = pltutils.get1Dvar(simname, "canopy", "rsshd")
   rsshd     = arsshd[:, tslice]

   # ansun
   z, aansun = pltutils.get1Dvar(simname, "canopy", "anetsun")
   ansun     = aansun[:, tslice] 

   # anshd
   z, aanshd = pltutils.get1Dvar(simname, "canopy", "anetshd")
   anshd     = aanshd[:, tslice] 

   # hsun
   z, ahsun = pltutils.get1Dvar(simname, "canopy", "hsun")
   hsun     = ahsun[:, tslice] 

   # hshd
   z, ahshd = pltutils.get1Dvar(simname, "canopy", "hshd")
   hshd     = ahshd[:, tslice] 

   # htot
   z, ahtot = pltutils.get1Dvar(simname, "canopy", "htot")
   htot     = ahtot[:, tslice] 

   # esun
   z, aesun = pltutils.get1Dvar(simname, "canopy", "esun")
   esun     = aesun[:, tslice] 

   # eshd
   z, aeshd = pltutils.get1Dvar(simname, "canopy", "eshd")
   eshd     = aeshd[:, tslice] 

   # etot
   z, aetot = pltutils.get1Dvar(simname, "canopy", "etot")
   etot     = aetot[:, tslice] 

   # gb
   z, agb   = pltutils.get1Dvar(simname, "canopy", "gb")
   gb       = agb[:, tslice] 

   # gl_sun
   z, agl_sun = pltutils.get1Dvar(simname, "canopy", "gl_sun")
   gl_sun     = agl_sun[:, tslice] 

   # gl_shd
   z, agl_shd = pltutils.get1Dvar(simname, "canopy", "gl_shd")
   gl_shd     = agl_shd[:, tslice] 

   # vl_sun
   z, avl_sun = pltutils.get1Dvar(simname, "canopy", "vl_sun")
   vl_sun     = avl_sun[:, tslice] 

   # vl_shd
   z, avl_shd = pltutils.get1Dvar(simname, "canopy", "vl_shd")
   vl_shd     = avl_shd[:, tslice] 


   # create the plots
   fig = plt.figure(figsize=(18,12))

   # sun/shade fractions
   ax1 = fig.add_subplot(2,6,1)
   pfsun, = ax1.plot(fsun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="f$_{sun}$")
   pfshd, = ax1.plot(fshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="f$_{shd}$")
   pltutils.setstdfmts(ax1, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.85,0.05))
   plt.xlim(xmax=1.0, xmin=0.0)
   plt.ylim(ymax=35.0, ymin=0.0)

   # Rabs sun/shade
   ax2 = fig.add_subplot(2,6,2)
   prasun, = ax2.plot(rasun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="R$_{a,sun}$")
   prashd, = ax2.plot(rashd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="R$_{a,shd}$")
   pltutils.setstdfmts(ax2, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.55,0.05))
   plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=35.0, ymin=0.0)
   plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   ax2.axes.yaxis.set_ticklabels([])

   # LW up & dn
   ax3 = fig.add_subplot(2,6,3)
   plwup, = ax3.plot(lwup, z, color=colors[4], linestyle="-", marker="None", linewidth=lnwdth, label="LW$_{up}$")
   plwdn, = ax3.plot(lwdn, z, color=colors[2], linestyle="-", marker="None", linewidth=lnwdth, label="LW$_{dn}$")
   pltutils.setstdfmts(ax3, tlmaj, tlmin, tlbsize, tlbpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
   plt.xlim(xmax=500.0, xmin=-20.0)
   plt.ylim(ymax=35.0, ymin=0.0)
   plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   ax3.axes.yaxis.set_ticklabels([])

   # Anet, Net photosynthetic assimilation rate, sun & shade
   ax4 = fig.add_subplot(2,6,4)
   pansun, = ax4.plot(ansun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="A$_{n,sun}$")
   panshd, = ax4.plot(anshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="A$_{n,shd}$")
   pltutils.setstdfmts(ax4, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("$\mu$mol m$^{-2}$ s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
   plt.xlim(xmax=20., xmin=-1.)
   plt.ylim(ymax=35.0, ymin=0.0)
   ax4.axes.yaxis.set_ticklabels([])

   # rt ???  TODO
   ax5 = fig.add_subplot(2,6,5)
   prtsun, = ax5.plot(rtsun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="R$_{t,sun}$")
   prtshd, = ax5.plot(rtshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="R$_{t,shd}$")
   pltutils.setstdfmts(ax5, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
   plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=35.0, ymin=0.0)
   ax5.axes.yaxis.set_ticklabels([])

   # PPFD
   ax6 = fig.add_subplot(2,6,6)
   pppfdsun, = ax6.plot(ppfdsun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="PPFD$_{sun}$")
   pppfdshd, = ax6.plot(ppfdshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="PPFD$_{shd}$")
   pltutils.setstdfmts(ax6, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("$\mu$mol m$^{-2}$ s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
   plt.xlim(xmax=1200.0, xmin=-100.0)
   plt.ylim(ymax=35.0, ymin=0.0)
   ax6.axes.yaxis.set_ticklabels([])

   # NIR
   ax7 = fig.add_subplot(2,6,7)
   pnirsun, = ax7.plot(nirsun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="NIR$_{sun}$")
   pnirshd, = ax7.plot(nirshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="NIR$_{shd}$")
   pltutils.setstdfmts(ax7, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("W/m$^2$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.65,0.05))
   plt.xlim(xmax=500.0, xmin=-20.0)
   plt.ylim(ymax=35.0, ymin=0.0)

   # rs, Stomatal resistance, sun & shade
   ax8 = fig.add_subplot(2,6,8)
   prssun, = ax8.plot(rssun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="r$_{s,sun}$")
   prsshd, = ax8.plot(rsshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="r$_{s,shd}$")
#  xmin, xmax = ax8.get_xlim()
#  plt.xticks(np.arange(xmin, xmax+1, 1000.))
   pltutils.setstdfmts(ax8, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("s m$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
   plt.xlim(xmax=2000., xmin=-100.)
   plt.ylim(ymax=35.0, ymin=0.0)
   ax8.axes.yaxis.set_ticklabels([])

   # canopy sensible heat flux, sun & shade
   ax9 = fig.add_subplot(2,6,9)
   phsun, = ax9.plot(hsun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="hsun")
   phshd, = ax9.plot(hshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="hshd")
   phtot, = ax9.plot(htot, z, color=colors[2], linestyle="-", marker="None", linewidth=lnwdth, label="htot")
   pltutils.setstdfmts(ax9, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("W cm$^{-2}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
#  plt.xlim(xmax=2000., xmin=-100.)
   plt.ylim(ymax=35.0, ymin=0.0)
   ax9.axes.yaxis.set_ticklabels([])

   # canopy evapotranspiration flux, sun & shade
   axa = fig.add_subplot(2,6,10)
   pesun, = axa.plot(esun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="esun")
   peshd, = axa.plot(eshd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="eshd")
   petot, = axa.plot(etot, z, color=colors[2], linestyle="-", marker="None", linewidth=lnwdth, label="etot")
   pltutils.setstdfmts(axa, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("mol cm$^{-2}$ s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
#  plt.xlim(xmax=2000., xmin=-100.)
   plt.ylim(ymax=35.0, ymin=0.0)
   axa.axes.yaxis.set_ticklabels([])

   # leaf conductance, sun & shade, and leaf boundary layer conductance
   axb = fig.add_subplot(2,6,11)
   pglsun, = axb.plot(gl_sun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="gl_sun")
   pglshd, = axb.plot(gl_shd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="gl_shd")
   pgb, = axb.plot(gb, z, color=colors[4], linestyle="-", marker="None", linewidth=lnwdth, label="gb")
   pltutils.setstdfmts(axb, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("mol m$^{-2}$ s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
   plt.ylim(ymax=35.0, ymin=0.0)
   axb.axes.yaxis.set_ticklabels([])

   # leaf exchange velocity, sun & shade
   axc = fig.add_subplot(2,6,12)
   pvlsun, = axc.plot(vl_sun, z, color=colors[3], linestyle="-", marker="None", linewidth=lnwdth, label="vl_sun")
   pvlshd, = axc.plot(vl_shd, z, color=colors[0], linestyle="-", marker="None", linewidth=lnwdth, label="vl_shd")
   pltutils.setstdfmts(axc, tlmaj, tlmin, tlbsize, tlbpad)
   plt.xlabel("cm s$^{-1}$", fontsize=xfsize, labelpad=xtpad)
   plt.legend(loc=4, fontsize=lfsize, bbox_to_anchor=(0.95,0.05))
#  plt.xlim(xmax=2000., xmin=-100.)
   plt.ylim(ymax=35.0, ymin=0.0)
   axc.axes.yaxis.set_ticklabels([])

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

       pfsun.set_xdata(afsun[:, tslice])
       pfshd.set_xdata(afshd[:, tslice])
       xmin=0.0
       xmax=1.0
       ax1.set_xlim(xmax=xmax, xmin=xmin)

#      xmaxsun = np.amax(arasun[:, tslice])
#      xminsun = np.amin(arasun[:, tslice])
#      xmaxshd = np.amax(arashd[:, tslice])
#      xminshd = np.amin(arashd[:, tslice])       
#      xmax = max(xmaxsun, xmaxshd)
#      xmin = min(xminsun, xminshd)
       xmax = 1200.0
       xmin = -100.0
       prasun.set_xdata(arasun[:, tslice])
       prashd.set_xdata(arashd[:, tslice])
       ax2.set_xlim(xmax=xmax, xmin=xmin)
#      ax2.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

       xmaxup = np.amax(alwup[:, tslice])
       xminup = np.amin(alwup[:, tslice])
       xmaxdn = np.amax(alwdn[:, tslice])
       xmindn = np.amin(alwdn[:, tslice])
       xmax = max(xmaxup, xmaxdn)
       xmin = min(xminup, xmindn)
       plwup.set_xdata(alwup[:, tslice])
       plwdn.set_xdata(alwdn[:, tslice])
       ax3.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

#      xmaxsun = np.amax(aansun[:, tslice])
#      xminsun = np.amin(aansun[:, tslice])
#      xmaxshd = np.amax(aanshd[:, tslice])
#      xminshd = np.amin(aanshd[:, tslice])
#      xmax = max(xmaxsun, xmaxshd)
#      xmin = min(xminsun, xminshd)
       xmax = 20.0
       xmin = -1.0
       pansun.set_xdata(aansun[:, tslice])
       panshd.set_xdata(aanshd[:, tslice])
#      ax4.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)
       ax4.set_xlim(xmax=xmax, xmin=xmin)

#      xmaxsun = np.amax(artsun[:, tslice])
#      xminsun = np.amin(artsun[:, tslice])
#      xmaxshd = np.amax(artshd[:, tslice])
#      xminshd = np.amin(artshd[:, tslice])
#      xmax = max(xmaxsun, xmaxshd)
#      xmin = min(xminsun, xminshd)
       xmax = 500.0
       xmin = -20.0
       prtsun.set_xdata(artsun[:, tslice])
       prtshd.set_xdata(artshd[:, tslice])
#      ax5.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)
       ax5.set_xlim(xmax=xmax, xmin=xmin)

       xmax = 2500.
       xmin = -20.0
       pppfdsun.set_xdata(appfdsun[:, tslice])
       pppfdshd.set_xdata(appfdshd[:, tslice])
       ax6.set_xlim(xmax=xmax, xmin=xmin)

#      xmaxsun = np.amax(anirsun[:, tslice])
#      xminsun = np.amin(anirsun[:, tslice])
#      xmaxshd = np.amax(anirshd[:, tslice])
#      xminshd = np.amin(anirshd[:, tslice])
#      xmax = max(xmaxsun, xmaxshd)
#      xmin = min(xminsun, xminshd)
       xmax = 500.0
       xmin = -20.0
       pnirsun.set_xdata(anirsun[:, tslice])
       pnirshd.set_xdata(anirshd[:, tslice])
#      ax7.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)
       ax7.set_xlim(xmax=xmax, xmin=xmin)

#      xmaxsun = np.amax(arssun[:, tslice])
#      xminsun = np.amin(arssun[:, tslice])
#      xmaxshd = np.amax(arsshd[:, tslice])
#      xminshd = np.amin(arsshd[:, tslice])
#      xmax = max(xmaxsun, xmaxshd)
#      xmin = min(xminsun, xminshd)
       xmax = 2000.0
       xmin = -100.0
       prssun.set_xdata(arssun[:, tslice])
       prsshd.set_xdata(arsshd[:, tslice])
#      ax8.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)
       ax8.set_xlim(xmax=xmax, xmin=xmin)

       xmax = np.amax(ahtot[:, tslice])
       xmin = np.amin(ahtot[:, tslice])
       phsun.set_xdata(ahsun[:, tslice])
       phshd.set_xdata(ahshd[:, tslice])
       phtot.set_xdata(ahtot[:, tslice])
       ax9.set_xlim(xmax=xmax*1.1, xmin=xmin)

       xmax = np.amax(aetot[:, tslice])
       xmin = np.amin(aetot[:, tslice])
       pesun.set_xdata(aesun[:, tslice])
       peshd.set_xdata(aeshd[:, tslice])
       petot.set_xdata(aetot[:, tslice])
       axa.set_xlim(xmax=xmax*1.1, xmin=xmin)

       xmaxsun = np.amax(agl_sun[:, tslice])
       xminsun = np.amin(agl_sun[:, tslice])
       xmaxshd = np.amax(agl_shd[:, tslice])
       xminshd = np.amin(agl_shd[:, tslice])
       xmaxgb = np.amax(agb[:, tslice])
       xmingb = np.amin(agb[:, tslice])
       xmaxgl = max(xmaxsun, xmaxshd)
       xmingl = min(xminsun, xminshd)
       xmax = max(xmaxgl, xmaxgb)
       xmin = min(xmingl, xmingb)
       pglsun.set_xdata(agl_sun[:, tslice])
       pglshd.set_xdata(agl_shd[:, tslice])
       pgb.set_xdata(agb[:, tslice])
       axb.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

       xmaxsun = np.amax(avl_sun[:, tslice])
       xminsun = np.amin(avl_sun[:, tslice])
       xmaxshd = np.amax(avl_shd[:, tslice])
       xminshd = np.amin(avl_shd[:, tslice])
       xmax = max(xmaxsun, xmaxshd)
       xmin = min(xminsun, xminshd)
       pvlsun.set_xdata(avl_sun[:, tslice])
       pvlshd.set_xdata(avl_shd[:, tslice])
       axc.set_xlim(xmax=xmax*1.1, xmin=xmin*0.95)

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
