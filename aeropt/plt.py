###########################################################################################
# aeropt/plt.py
#
# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
#
# Author:
#    Ramiro Checa-Garcia. ECMWF
#
# Modifications:
#    10-Jan-2023   Ramiro Checa-Garcia    1st. version
#
#                                                                                         
# Info: 
#      Provides CLASSES and FUNCTIONS:
#                                                                                         
#      FUNCTIONS        
#       * aerplt               :
#       * plt_csv              :
#       * plt_csv_refind       :
#       * plt_nc               :
#       * plt_nc_refindex      :
#       * one_col              :
#       * few_col              :
#       * phase_col            :
#       * phase_one            :
#       * plt_nc_optical       : 
#       * plt_nc_phasefunction :
##########################################################################################



import xarray as xr
import numpy as np

try:
    import tomllib as toml
except ModuleNotFoundError:
    import toml

import matplotlib.pyplot as plt
import os
import glob 

plt.rcParams.update({'axes.labelsize': '9'})

def aerplt(output, data, extra=dict()):
    if "pltdir" not in extra.keys():
       extra["pltdir"]="outputplt"
    
    if os.path.isfile(data):
        if data[-3::]==".nc":
            plt_nc(output, data, extra=extra)
        else:
            plt_csv(output, data, extra=extra)
    return

def plt_csv(output, data, extra):

    if output=="refindex":
        plt_csv_refindex(data, extra)

    return

def plt_csv_refind(data, extra):

    return


def plt_nc(output, data, extra):

    if output=="refindex":
        plt_nc_refindex(data, extra)
    elif output=='optical':
        plt_nc_optical(data, extra)
    elif output=="phasefunction":
        plt_nc_phasefunction(data, extra)
    return

def plt_nc_refindex(data, extra, pformat="png"):
    ds = xr.open_dataset(data)
    pname = "_".join(os.path.basename(data).split('_')[0:-1])
    sname = os.path.join(extra["pltdir"],"refindex_"+pname+'.'+pformat)
    if "source_config" in ds.attrs.keys():
        varnames=["ref_idx_real","ref_idx_img"]
        fig, axs = plt.subplots(nrows=2, sharex=True)
        fig.set_size_inches(5.5, 6.5)
        fig.suptitle(pname, x=0.1, y=0.985, ha='left', fontsize=10)
        one_col(axs, ds, varnames)

        plt.subplots_adjust(left=0.2)
        #plt.tight_layout()
        plt.savefig(sname)
        plt.close(fig)
    else:
        varnames=[]

    return

def one_col(axs, ds, varnames):
    for ii in range(len(axs)):
        if ii!=0:
            addlegend=False
        else:
            addlegend=True
        if 'size_bin' in ds[varnames[ii]].dims:
           ds[varnames[ii]].isel(size_bin=0).plot.line(ax=axs[ii], xscale='log', 
                                                       hue="rel_hum", add_legend=addlegend, **{'lw':0.5})
        else:
           ds[varnames[ii]].plot.line(ax=axs[ii], xscale='log', hue="rel_hum", 
                                      add_legend=addlegend, **{'lw':0.5})
        axs[ii].set_title("")
        if ii != len(axs)-1:
                axs[ii].set_xlabel("")
        if ii == 0:
                handles, labels = axs[ii].get_legend_handles_labels()
                labels=ds["rel_hum"][:].values
                axs[ii].legend(labels, ncol=6, fontsize=7, title='rel_hum[%]', 
                            title_fontsize=8 , loc='lower right',
                            frameon=False, bbox_to_anchor=(1.0, 1.01))
    return 

def few_col(axs, ds, varnames):
    axs0=axs.shape[0]
    axs1=axs.shape[1]
    for ii in range(axs0):
        for jj in range(axs1):
            if ii==0 and jj==axs1-1:
                addlegend=True
            else:
                addlegend=False
            ds[varnames[ii]].isel(size_bin=jj).plot.line(ax=axs[ii,jj], xscale='log', 
                                                         hue="rel_hum", add_legend=addlegend, **{'lw':0.5})
            axs[ii,jj].set_title("")
            if ii!=axs0-1:
                axs[ii,jj].set_xlabel("")
            if jj!=0:
                axs[ii,jj].set_ylabel("")
            if ii == 0 and jj==axs1-1:
                handles, labels = axs[ii,jj].get_legend_handles_labels()
                labels=ds["rel_hum"][:].values
                axs[ii,jj].legend(labels, ncol=6, fontsize=7, title='rel_hum[%]', 
                                    title_fontsize=8 , loc='lower right',
                                    frameon=False, bbox_to_anchor=(1.0, 1.01))

    return

def phase_col(axs, ds, varname):
    axs0=axs.shape[0]
    axs1=axs.shape[1]
    wl  = 1.e9*ds["wavelength"][:].values
    iwl = [0,8,14]
    for ii in range(axs0):
        for jj in range(axs1):
            if ii==0 and jj==axs1-1:
                addlegend=True
            else:
                addlegend=False
            dtplot = ds[varname[0]].isel(size_bin=ii, wavelength=iwl[jj]).plot.line(ax=axs[ii,jj], 
                                                                                    yscale='log', 
                                                                                    hue="rel_hum",
                                                                                    add_legend=addlegend, **{'lw':0.5})
            stitle = " bin ="+str(ii)+"\n wl ="+str(wl[iwl[jj]])+' nm'
            axs[ii,jj].text(0.99, 0.95, stitle, horizontalalignment='right', 
                            verticalalignment='top', transform=axs[ii,jj].transAxes)
            axs[ii,jj].set_title("")
            if ii!=axs0-1:
                axs[ii,jj].set_xlabel("")
            if jj!=0:
                axs[ii,jj].set_ylabel("")
            if ii == 0 and jj==axs1-1:
                handles, labels = axs[ii,jj].get_legend_handles_labels()
                labels=ds["rel_hum"][:].values
                axs[ii,jj].legend(labels, ncol=6, fontsize=7, title='rel_hum[%]', 
                                    title_fontsize=8 , loc='lower right',
                                    frameon=False, bbox_to_anchor=(1.0, 1.01))
    return

def phase_one(axs, ds, varname):
    axs0=axs.shape[0]
    wl  = 1.e9*ds["wavelength"][:].values
    iwl = [0,8,14]
    ii  = 0
    for jj in range(axs0):
            if ii==0 and jj==axs0-1:
                addlegend=True
            else:
                addlegend=False
            dtplot = ds[varname[0]].isel(size_bin=ii, wavelength=iwl[jj]).plot.line(ax=axs[jj],
                                                                                    yscale='log',
                                                                                    hue="rel_hum",
                                                                                    add_legend=addlegend, **{'lw':0.5})
            stitle = str(wl[iwl[jj]])+' nm'
            axs[jj].text(0.99, 0.97, stitle, horizontalalignment='right', verticalalignment='top', transform=axs[jj].transAxes)
            axs[jj].set_title("")
            if jj!=0:
                axs[jj].set_ylabel("")
            if ii == 0 and jj==axs0-1:
                handles, labels = axs[jj].get_legend_handles_labels()
                labels=ds["rel_hum"][:].values
                axs[jj].legend(labels, ncol=6, fontsize=7, title='rel_hum[%]', 
                                    title_fontsize=8 , loc='lower right',
                                    frameon=False, bbox_to_anchor=(1.0, 1.01))
    return




def plt_nc_optical(data, extra, pformat="png"):
    ds  = xr.open_dataset(data)
    pname = "_".join(os.path.basename(data).split('_')[0:-1])
    sname = os.path.join(extra["pltdir"],"optproperties_"+pname+'.'+pformat)
    if "source_config" in ds.attrs.keys():
        ncol=len(ds['size_bin'][:].values)
        varnames=["extinction","single_scatter_albedo","asymmetry_factor"]
        fig, axs = plt.subplots(nrows=3, ncols=ncol, sharex=True, sharey='row')
        fig.set_size_inches(5+(ncol-1)*2.5, 8)
        fig.suptitle(pname, x=0.1, y=0.985, ha='left', fontsize=10)
        if ncol==1:
            one_col(axs, ds, varnames)
        if ncol>1:
            few_col(axs, ds, varnames)
        plt.subplots_adjust(left=0.2)
        #plt.tight_layout()
        plt.savefig(sname)
        plt.close(fig)
    else:
        varnames=[]
    return

def plt_nc_phasefunction(data, extra, pformat="png"):
    ds  = xr.open_dataset(data)
    pname = "_".join(os.path.basename(data).split('_')[0:-1])
    sname = os.path.join(extra["pltdir"],"phasefunction_"+pname+'.'+pformat)
    if "source_config" in ds.attrs.keys():
        ncol=len(ds['size_bin'][:].values)
        varnames=["phase_function"]
        if 'size_bin' in ds[varnames[0]].dims:
            if ncol > 1:
               fig, axs = plt.subplots(nrows=3, ncols=ncol, sharex=True, sharey='row')
               fig.set_size_inches(5+(ncol-1)*2.5, 8)
            if ncol == 1:
               fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True, sharey='row')
               fig.set_size_inches(5+(3-1)*2.5, 3.5)

        fig.suptitle(pname, x=0.1, y=0.985, ha='left', fontsize=10)
        if ncol > 1:
            phase_col(axs, ds, varnames)
            plt.subplots_adjust(left=0.2)
        else:
            phase_one(axs, ds, varnames)
            plt.subplots_adjust(left=0.2, top=0.75, bottom=0.2)
        plt.savefig(sname)
        plt.close(fig)
    else:
        varnames=[]
    return



