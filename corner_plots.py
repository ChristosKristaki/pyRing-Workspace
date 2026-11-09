import numpy as np
import corner
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import glob as glob
import re
import h5py
from numpy.lib import recfunctions as rfn
import configparser 
from matplotlib.colors import to_rgba
import os
import copy

STYLE = dict(
    bins=50,
    smooth=1.5,
    smooth1d=1.0,
    levels=[0.5, 0.9],
    quantiles=[0.5],
    true_quantiles=[0.05, 0.5, 0.95],
    title_kwargs={'fontsize': 12},
    use_math_text=True,
    max_n_ticks=3,
    plot_datapoints=False,
    plot_density=False,
    plot_contours=True,
    fill_contours=True,
    no_fill_contours=True,
    contour_kwargs={'linewidths': 1.5},
    hist_kwargs={'linewidth': 1.5},
)

def ramp(color, top=0.3):
    n = len(STYLE['levels']) + 1
    return [to_rgba(color, top * i / (n - 1)) for i in range(n)]

def suffix_key(f):
    m = re.search(r'rin_S\d{6}([a-z]*)_', f)
    s = m.group(1) if m else ''
    return (len(s), s) 

def w(d):
    return np.ones(len(d)) / len(d)


def get_median_and_errors(data, labels):
    results = {}
    for i, label in enumerate(labels):
        lower, median, upper = np.percentile(data[:, i], [5, 50, 95])
        results[label] = (median, median - lower, upper - median)
    return results

def get_median_and_errors_single(data):
    lower, median, upper = np.percentile(data, [5, 50, 95])
    results = (median, median - lower, upper - median)
    return results

def corner_plots_220(gw_event, dir_aq_220, dir_gr_220, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_220, skiprows=1)
    data_gr = np.loadtxt(dir_gr_220, skiprows=1)

    data_aq_shared = np.delete(data_aq, 4, axis=1)
    data_aq_physical = data_aq_shared[:, :-2]
    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_aq = get_median_and_errors(data_aq_physical, labels)
    results_gr = get_median_and_errors(data_gr_physical, labels)

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue', truths = truths,
                        truth_color='k', contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE))
    corner.corner(data_aq_physical, labels=labels, color='red', fig=figure, 
                        contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GRvsAQ 220 mode {start_time} ', fontsize=20, y=0.7, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_aq[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"AQ: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    results_snr_gr_220, results_snr_aq_220 = get_snrs(gw_event, start_time, gr220=True, aq220=True, gr221=False, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label='GR run with high settings and 220 mode'),
            mlines.Line2D([], [], color='red', label='AQ run with high settings and 220 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_220[0]:.2f}$_{{-{results_snr_gr_220[1]:.2f}}}^{{+{results_snr_gr_220[2]:.2f}}}$'),
            mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_220[0]:.2f}$_{{-{results_snr_aq_220[1]:.2f}}}^{{+{results_snr_aq_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)

    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/AQvsGR_220/220_GRvsAQ_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/AQvsGR_220/220_GRvsAQ_{gw_event}.png', dpi=150)

    plt.close(figure)

def corner_plots_221(gw_event, dir_aq_221, dir_gr_221, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_221, skiprows=1)
    data_gr = np.loadtxt(dir_gr_221, skiprows=1)

    data_aq_shared = np.delete(data_aq, 4, axis=1)
    data_aq_physical = data_aq_shared[:, :-2]
    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'A2221', 'phi2220', 'phi2221']

    results_aq = get_median_and_errors(data_aq_physical, labels)
    results_gr = get_median_and_errors(data_gr_physical, labels)

    truths = [None, mf_truth, None, None, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue', truths = truths,
                            truth_color='k',  contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE))
    corner.corner(data_aq_physical, labels=labels, color='red', fig=figure, 
                             contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GRvsAQ 221 mode {start_time}', fontsize=30, y=0.7, x=0.7, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_aq[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"AQ: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    results_snr_gr_221, results_snr_aq_221 = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=True, aq221=True)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label='GR run with high settings and 221 mode'),
            mlines.Line2D([], [], color='red', label='AQ run with high settings and 221 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_221[0]:.2f}$_{{-{results_snr_gr_221[1]:.2f}}}^{{+{results_snr_gr_221[2]:.2f}}}$'),
            mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_221[0]:.2f}$_{{-{results_snr_aq_221[1]:.2f}}}^{{+{results_snr_aq_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=25)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=25)
   

    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/AQvsGR_221/220+221_GRvsAQ_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/AQvsGR_221/220+221_GRvsAQ_{gw_event}.png', dpi=150)

    plt.close(figure)


def corner_plots_GWTC4_220(gw_event, dir_aq_220, dir_gr_220, dir_gwtc, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_220, skiprows=1)
    data_gr = np.loadtxt(dir_gr_220, skiprows=1)

    if gw_catalogue_flag == 'GWTC3':
        with h5py.File(dir_gwtc[0], "r") as f:
            datah5_gwtc = f['GWTC3_rerun_PROD1/posterior_samples'][()]
            data_gwtc = rfn.structured_to_unstructured(datah5_gwtc)          # 2D: (n_samples, n_fields)
            data_gwtc_shared = np.delete(data_gwtc, 5, axis=1)
            data_gwtc_physical = data_gwtc_shared[:, :-7]
    elif gw_catalogue_flag == 'GWTC4':
        data_gwtc = np.loadtxt(dir_gwtc, skiprows=1)
        data_gwtc_physical = data_gwtc[:, :-1]


    data_aq_shared = np.delete(data_aq, 4, axis=1)
    data_aq_physical = data_aq_shared[:, :-2]
    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    #print(dir_gr_220)
    #print(os.path.getmtime(dir_gr_220), data_gr_physical.shape)
    #print(get_median_and_errors(data_gr_physical[:, [1, 2]], ['Mf', 'af']))

    results_aq = get_median_and_errors(data_aq_physical, labels)
    results_gr = get_median_and_errors(data_gr_physical, labels)
    results_gwtc = get_median_and_errors(data_gwtc_physical, labels)

    # Compute shared ranges across both datasets so bins line up
    combined = np.vstack([data_gr_physical, data_gwtc_physical, data_aq_physical])
    ranges = [(combined[:, i].min(), combined[:, i].max()) for i in range(combined.shape[1])]

    print(len(combined))
    print(ranges)

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical,labels=labels, color='blue', 
                        truths = truths,
                        contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE), weights=w(data_gr_physical)
                        )
    corner.corner(data_gwtc_physical, labels=labels, color='green',
                        fig = figure,
                        contourf_kwargs={'colors': ramp('green')}, **copy.deepcopy(STYLE), weights=w(data_gwtc_physical)
                        )
    corner.corner(data_aq_physical, labels=labels, color='red',
                        fig = figure,
                        contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE), weights=w(data_aq_physical))

    #normalization check
    for i in data_gr_physical, data_aq_physical, data_gwtc_physical:
        h, edges = np.histogram(i[:, 1], bins=50, range=ranges[1], density=True)
        print(np.sum(h * np.diff(edges)))     # should print 1.0

    for name, d in [("GR", data_gr_physical), ("AQ", data_aq_physical), ("GWTC", data_gwtc_physical)]: print(name, d.shape, "NaNs:", int(np.isnan(d).any(axis=1).sum()))

    #figure = corner.corner(data_gr_physical, labels=labels, color='blue',
    #                     truth_color='blue')
    #corner.corner(data_aq_physical, labels=labels, color='red', fig=figure,
    #           truth_color='red')
    #corner.corner(data_gwtc_physical, labels=labels, color='green', fig=figure,
    #           truth_color='green')

    figure.suptitle(f'GW{gw_event} GRvsAQvsGWTC4 220 mode {start_time}', fontsize=20, y=0.75, x=0.7, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_aq[label]
        med3, minus3, plus3 = results_gwtc[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"AQ: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$\n"
                f"GWTC4: {med3:.2f}$_{{-{minus3:.2f}}}^{{+{plus3:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    results_snr_gr_220, results_snr_aq_220 = get_snrs(gw_event, start_time, gr220=True, aq220=True, gr221=False, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label=f'GR run with high settings and 220 mode at {start_time}'),
            mlines.Line2D([], [], color='red', label=f'AQ run with high settings and 220 mode at {start_time}'),
            mlines.Line2D([], [], color='green', label='GWTC samples for 220 mode at 10M')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_220[0]:.2f}$_{{-{results_snr_gr_220[1]:.2f}}}^{{+{results_snr_gr_220[2]:.2f}}}$'),
            mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_220[0]:.2f}$_{{-{results_snr_aq_220[1]:.2f}}}^{{+{results_snr_aq_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)

    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/AQvsGRvsGWTC3/220_GRvsAQvsGWTC3_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/AQvsGRvsGWTC4/220_GRvsAQvsGWTC4_{gw_event}.png', dpi=150)
    
    plt.close(figure)


def corner_plots_GWTC4_221(gw_event, dir_aq_221, dir_gr_221, dir_gwtc, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_221, skiprows=1)
    data_gr = np.loadtxt(dir_gr_221, skiprows=1)
    
    if gw_catalogue_flag == 'GWTC3':
        with h5py.File(dir_gwtc[0], "r") as f:
            datah5_gwtc = f['GWTC3_rerun_PROD1/posterior_samples'][()]
            data_gwtc = rfn.structured_to_unstructured(datah5_gwtc)          # 2D: (n_samples, n_fields)
            data_gwtc_shared = np.delete(data_gwtc, 5, axis=1)
            data_gwtc_physical = data_gwtc_shared[:, :-7]
    elif gw_catalogue_flag == 'GWTC4':
        data_gwtc = np.loadtxt(dir_gwtc, skiprows=1)
        data_gwtc_physical = data_gwtc[:, :-1]

    data_aq_shared = np.delete(data_aq, [4,6,8], axis=1)
    data_gr_shared = np.delete(data_gr, [5,7], axis=1)
    data_aq_physical = data_aq_shared[:, :-2]
    data_gr_physical = data_gr_shared[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_aq = get_median_and_errors(data_aq_physical, labels)
    results_gr = get_median_and_errors(data_gr_physical, labels)
    results_gwtc = get_median_and_errors(data_gwtc_physical, labels)

    # Compute shared ranges across both datasets so bins line up
    combined = np.vstack([data_gr_physical, data_gwtc_physical, data_aq_physical])
    ranges = [(combined[:, i].min(), combined[:, i].max()) for i in range(combined.shape[1])]
    
    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical,labels=labels, color='blue',
                        contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE),
                        truths = truths,weights=w(data_gr_physical)
                        )
    corner.corner(data_gwtc_physical, labels=labels, color='green',
                        fig = figure,weights=w(data_gwtc_physical),
                        contourf_kwargs={'colors': ramp('green')}, **copy.deepcopy(STYLE)
                        )
    corner.corner(data_aq_physical, labels=labels, color='red',weights=w(data_aq_physical),
                        fig = figure, contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE)
                        )

    #normalization check
    for i in data_gr_physical, data_aq_physical, data_gwtc_physical:
        h, edges = np.histogram(i[:, 1], bins=30, range=ranges[1], density=True)
        print(np.sum(h * np.diff(edges)))     # should print 1.0

    for name, d in [("GR", data_gr_physical), ("AQ", data_aq_physical), ("GWTC", data_gwtc_physical)]: print(name, d.shape, "NaNs:", int(np.isnan(d).any(axis=1).sum()))

    figure.suptitle(f'GW{gw_event} GRvsAQvsGWTC4 220 mode {start_time}', fontsize=20, y=0.75, x=0.7, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_aq[label]
        med3, minus3, plus3 = results_gwtc[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"AQ: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$\n"
                f"GWTC4: {med3:.2f}$_{{-{minus3:.2f}}}^{{+{plus3:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    results_snr_gr_221, results_snr_aq_221 = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=True, aq221=True)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label=f'GR run with high settings and 220+221 mode at {start_time}'),
            mlines.Line2D([], [], color='red', label=f'AQ run with high settings and 220+221 mode at {start_time}'),
            mlines.Line2D([], [], color='green', label='GWTC samples for 220 mode at 10M')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_221[0]:.2f}$_{{-{results_snr_gr_221[1]:.2f}}}^{{+{results_snr_gr_221[2]:.2f}}}$'),
            mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_221[0]:.2f}$_{{-{results_snr_aq_221[1]:.2f}}}^{{+{results_snr_aq_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)

    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/AQvsGRvsGWTC3/221_GRvsAQvsGWTC3_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/AQvsGRvsGWTC4/221_GRvsAQvsGWTC4_{gw_event}.png', dpi=150)
    
    plt.close(figure)


def cp_mfvsalpha_220(gw_event, dir_aq_220, gw_catalogue_flag, mf_truth, start_time):

    data = np.loadtxt(dir_aq_220, skiprows=1)

    mf    = data[:, 1]   # Mf
    alpha = data[:, 4]   # alpha

    samples = np.column_stack([mf, alpha])
    labels  = ['Mf', 'alpha']

    results = get_median_and_errors(samples, labels)
    
    truths = [mf_truth, None]


    figure = corner.corner(samples, labels=labels, color='red',
                            show_titles=True, title_fmt='.2f',
                            plot_datapoints=True, plot_density=True, plot_contours=True,
                            truths = truths,
                            truth_color='k', levels=[0.9], quantiles=[0.05, 0.95], title_quantiles=[0.05, 0.5, 0.95])

    w, h = figure.get_size_inches()
    figure.set_size_inches(w, h * 1.15)   # grow the canvas to make room

    figure.suptitle(f'GW{gw_event} Mf vs alpha, 220 mode {start_time}', fontsize=15, y=0.98, color='black')
    figure.subplots_adjust(top=0.90)      # push panels down to leave gap under the title

    results_snr_aq_220, = get_snrs(gw_event, start_time, gr220=False, aq220=True, gr221=False, aq221=False)

    #Legends

    snr =  [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_220[0]:.2f}$_{{-{results_snr_aq_220[1]:.2f}}}^{{+{results_snr_aq_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 2: right edge aligned with legend 1, directly underneath it
    figure.legend(handles=snr, loc = 'upper right', bbox_to_anchor=(0.95, 0.9), fontsize=10) 
 
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/Mfvsalpha_220/{gw_event}_Mfvsalpha_220.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/Mfvsalpha_220/{gw_event}_Mfvsalpha_220.png', dpi=150)

    plt.close(figure)



def cp_mfvsalpha_221(gw_event, dir_aq_221, gw_catalogue_flag, mf_truth, start_time):

    data = np.loadtxt(dir_aq_221, skiprows=1)

    mf    = data[:, 1]   # Mf
    alpha = data[:, 4]   # alpha

    samples = np.column_stack([mf, alpha])
    labels  = ['Mf', 'alpha']

    results = get_median_and_errors(samples, labels)

    truths = [mf_truth, None]

    figure = corner.corner(samples, labels=labels, color='red',
                            show_titles=True, title_fmt='.2f',
                            plot_datapoints=True, plot_density=True, plot_contours=True,
                            truths = truths,
                            truth_color='k', levels=[0.9], quantiles=[0.05, 0.95], title_quantiles=[0.25, 0.5, 0.75],   # title shows median with the 90% bounds
 )

    w, h = figure.get_size_inches()
    figure.set_size_inches(w, h * 1.15)   # grow the canvas to make room

    figure.suptitle(f'GW{gw_event} Mf vs alpha, 220+221 mode {start_time}', fontsize=15, y=0.98, color='black')
    figure.subplots_adjust(top=0.90)      # push panels down to leave gap under the title
    
    results_snr_aq_221, = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=False, aq221=True)

    #Legends

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_221[0]:.2f}$_{{-{results_snr_aq_221[1]:.2f}}}^{{+{results_snr_aq_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 2: right edge aligned with legend 1, directly underneath it
    figure.legend(handles=snr, loc = 'upper right', bbox_to_anchor=(0.95, 0.9), fontsize=10) 

    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/Mfvsalpha_221/{gw_event}_Mfvsalpha_220+221.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/Mfvsalpha_221/{gw_event}_Mfvsalpha_220+221.png', dpi=150)
    
    plt.close(figure)


def corner_plots_GRvsGWTC_220(gw_event, dir_gr_220, dir_gwtc, gw_catalogue_flag, mf_truth, start_time):

    data_gr = np.loadtxt(dir_gr_220, skiprows=1)

    if gw_catalogue_flag == 'GWTC3':
        with h5py.File(dir_gwtc[0], "r") as f:
            datah5_gwtc = f['GWTC3_rerun_PROD1/posterior_samples'][()]
            data_gwtc = rfn.structured_to_unstructured(datah5_gwtc)          # 2D: (n_samples, n_fields)
            data_gwtc_shared = np.delete(data_gwtc, 5, axis=1)
            data_gwtc_physical = data_gwtc_shared[:, :-7]
    elif gw_catalogue_flag == 'GWTC4':
        data_gwtc = np.loadtxt(dir_gwtc, skiprows=1)
        data_gwtc_physical = data_gwtc[:, :-1]

    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_gwtc = get_median_and_errors(data_gwtc_physical, labels)
    results_gr = get_median_and_errors(data_gr_physical, labels)

    # Compute shared ranges across both datasets so bins line up
    combined = np.vstack([data_gr_physical, data_gwtc_physical])
    ranges = [(combined[:, i].min(), combined[:, i].max()) for i in range(combined.shape[1])]

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue', weights=w(data_gr_physical),
                        contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE), truths = truths,
    )
    corner.corner(data_gwtc_physical, labels=labels, color='green', fig=figure, weights=w(data_gwtc_physical),
                contourf_kwargs={'colors': ramp('green')}, **copy.deepcopy(STYLE))

    #figure = corner.corner(data_gr_physical, labels=labels, color='blue',
     #                    truth_color='blue')
    #corner.corner(data_gwtc_physical, labels=labels, color='green', fig=figure,
      #         truth_color='red')

    figure.suptitle(f'GW{gw_event} GRvsGWTC4 220 mode {start_time}', fontsize=25, y=0.8, x=0.7, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_gwtc[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"GWTC: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    results_snr_gr_220, = get_snrs(gw_event, start_time, gr220=True, aq220=False, gr221=False, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label=f'GR run with high settings and 220 mode at {start_time}'),
            mlines.Line2D([], [], color='green', label='GWTC samples for 220 mode at 10M')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_220[0]:.2f}$_{{-{results_snr_gr_220[1]:.2f}}}^{{+{results_snr_gr_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)

    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/GRvsGWTC3_220/220_GRvsGWTC3_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/GRvsGWTC4_220/220_GRvsGWTC4_{gw_event}.png', dpi=150)

    plt.close(figure)

def corner_plots_GRvsGWTC_221(gw_event, dir_gr_221, dir_gwtc, gw_catalogue_flag, mf_truth, start_time):

    data_gr = np.loadtxt(dir_gr_221, skiprows = 1)

    if gw_catalogue_flag == 'GWTC3':
        with h5py.File(dir_gwtc[0], "r") as f:
            datah5_gwtc = f['GWTC3_rerun_PROD1/posterior_samples'][()]
            data_gwtc = rfn.structured_to_unstructured(datah5_gwtc)          # 2D: (n_samples, n_fields)
            data_gwtc_shared = np.delete(data_gwtc, 5, axis=1)
            data_gwtc_physical = data_gwtc_shared[:, :-7]
    elif gw_catalogue_flag == 'GWTC4':
        data_gwtc = np.loadtxt(dir_gwtc, skiprows=1)
        data_gwtc_physical = data_gwtc[:, :-1]

    data_gr_physical = data_gr[:, :-2]

    temp_gr1 = np.delete(data_gr, [5,7], axis=1)

    data_gr_physical = temp_gr1[:, :-2]

    #print("Shape:", data_gr_physical.shape)
    #print("Shape:", data_gwtc_physical.shape)

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_gr = get_median_and_errors(data_gr_physical, labels)
    results_gwtc = get_median_and_errors(data_gwtc_physical, labels)

    #print("At corner call — data shape:", data_gr_physical.shape)
    #print("At corner call — labels:", labels, "len:", len(labels))

    # Compute shared ranges across both datasets so bins line up
    combined = np.vstack([data_gr_physical, data_gwtc_physical])
    ranges = [(combined[:, i].min(), combined[:, i].max()) for i in range(combined.shape[1])]

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue',
                        truths = truths, weights=w(data_gr_physical),
                        truth_color='k',  contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE))
    corner.corner(data_gwtc_physical, labels=labels, color='green', fig=figure, weights=w(data_gwtc_physical),
                        contourf_kwargs={'colors': ramp('green')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GR221 vs GWTC4 run {start_time}' ,fontsize=25, y=0.8, x=0.7, color='black')

     # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        med2, minus2, plus2 = results_gwtc[label]
        title = (f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"GWTC4: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')
    
    results_snr_gr_221, = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=True, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label=f'GR run with high settings and 220+221 mode at {start_time}'),
            mlines.Line2D([], [], color='green', label='GWTC samples for 220+221 mode at 10M')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_221[0]:.2f}$_{{-{results_snr_gr_221[1]:.2f}}}^{{+{results_snr_gr_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)
    
    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/GRvsGWTC3_221/221_GRvs221_GWTC3_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/GRvsGWTC4_221/221_GRvs221_GWTC4_{gw_event}.png', dpi=150)

    plt.close(figure)


def posterior_alpha_220(gw_event, dir_aq_220, gw_catalogue_flag, start_time):

    data_alpha = np.loadtxt(dir_aq_220, skiprows = 1)[:, 4]

    labels = {'alpha'}

    plt.hist(data_alpha, bins=50, color='red', edgecolor='black')

    results = get_median_and_errors_single(data_alpha)

    med, minus, plus = results

    plt.title(f'Posterior of the alpha parameter for the 220 only mode of {gw_event} \n a = {med:.2f}$_{{-{minus:.2f}}}^{{+{plus:.2f}}}$')
    plt.ylabel('Number of samples')
    plt.xlabel('Alpha')

    results_snr_aq_220, = get_snrs(gw_event, start_time, gr220=False, aq220=True, gr221=False, aq221=False)

    #Legends

    snr =  [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_220[0]:.2f}$_{{-{results_snr_aq_220[1]:.2f}}}^{{+{results_snr_aq_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    ax = plt.gca()

    # Legend 2: right edge aligned with legend 1, directly underneath it
    ax.legend(
        handles=snr,
        loc='upper left',
        bbox_to_anchor=(0, 1),
        fontsize=7
    )

    if gw_catalogue_flag == 'GWTC3':
        plt.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/alpha_posteriors/GW{gw_event}_alpha_220.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        plt.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/alpha_posteriors/GW{gw_event}_alpha_220.png', dpi=150)

    plt.close()

def posterior_alpha_221(gw_event, dir_aq_221, gw_catalogue_flag, start_time):

    data_alpha = np.loadtxt(dir_aq_221 , skiprows=1)[:,4]

    labels = {'alpha'}

    plt.hist(data_alpha, bins=50, color='red', edgecolor='black')

    results = get_median_and_errors_single(data_alpha)

    med, minus, plus = results

    plt.title(f'Posterior of the alpha parameter for the 220+221 mode of {gw_event} \n a = {med:.2f}$_{{-{minus:.2f}}}^{{+{plus:.2f}}}$')
    plt.ylabel('Number of samples')
    plt.xlabel('Alpha')


    ax = plt.gca()

    results_snr_aq_221, = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=False, aq221=True)

    #Legends

    snr =  [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_221[0]:.2f}$_{{-{results_snr_aq_221[1]:.2f}}}^{{+{results_snr_aq_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 2: right edge aligned with legend 1, directly underneath it
    ax.legend(
        handles=snr,
        loc='upper left',
        bbox_to_anchor=(0, 1),
        fontsize=7
    )

    if gw_catalogue_flag == 'GWTC3':
        plt.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/alpha_posteriors/GW{gw_event}_alpha_220_221.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        plt.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/alpha_posteriors/GW{gw_event}_alpha_220_221.png', dpi=150)

    plt.close()

def corner_plots_GWTC_only(gw_event, dir_gwtc, gw_catalogue_flag, mf_truth, start_time):

    if gw_catalogue_flag == 'GWTC3':
        with h5py.File(dir_gwtc[0], "r") as f:
            datah5_gwtc = f['GWTC3_rerun_PROD1/posterior_samples'][()]
            data_gwtc = rfn.structured_to_unstructured(datah5_gwtc)          # 2D: (n_samples, n_fields)
            data_gwtc_shared = np.delete(data_gwtc, 5, axis=1)
            data_gwtc_physical = data_gwtc_shared[:, :-7]
    elif gw_catalogue_flag == 'GWTC4':
        data_gwtc = np.loadtxt(dir_gwtc, skiprows=1)
        data_gwtc_physical = data_gwtc[:, :-1]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_gwtc = get_median_and_errors(data_gwtc_physical, labels)

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gwtc_physical, labels=labels, color='green',
                        truths = truths,
                        truth_color='k',  contourf_kwargs={'colors': ramp('green')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GWTC4 220 mode at 10M', fontsize=30, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))

    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med2, minus2, plus2 = results_gwtc[label]
        title = (f"GWTC: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')

    #Legend
    green_line = mlines.Line2D([], [], color='green', label='GWTC4 run with high settings and 220 mode')
    figure.legend(handles=[green_line], loc='upper right', fontsize=25)

    #saving the corner plot

    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/GWTC3_220/GWTC4_{gw_event}.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/GWTC4_220/GWTC4_{gw_event}.png', dpi=150)

    plt.close(figure)

def single_corner_plot_GR_220(gw_event, dir_gr_220, gw_catalogue_flag, mf_truth, start_time):

    data_gr = np.loadtxt(dir_gr_220, skiprows=1)

    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_gr = get_median_and_errors(data_gr_physical, labels)

    truths = [None, mf_truth, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue',
                        truths = truths,
                        truth_color='k',  contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GR 220 mode at {start_time}', fontsize=20, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        title = f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
        ax.set_title(title, fontsize=7, color='black')
    
    results_snr_gr_220, = get_snrs(gw_event, start_time, gr220=True, aq220=False, gr221=False, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label='GR run with high settings and 220 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_220[0]:.2f}$_{{-{results_snr_gr_220[1]:.2f}}}^{{+{results_snr_gr_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20) 
    
    #saving the corner plot

    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/full_posts/GW_{gw_event}_220_GR.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/full_posts/GW_{gw_event}_220_GR.png', dpi=150)
    
    plt.close(figure)

def single_corner_plot_AQ_220(gw_event, dir_aq_220, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_220, skiprows=1)

    data_aq_physical = data_aq[:, :-2]

    labels = ['psi', 'Mf', 'af', 'alpha', 'cosiota', 'A2220', 'phi2220']

    results_aq = get_median_and_errors(data_aq_physical, labels)

    truths = [None, mf_truth, None, None, None, None, None]

    figure = corner.corner(data_aq_physical, labels=labels, color='red',
                        truths = truths,
                        truth_color='k',  contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} AQ 220 mode at {start_time}', fontsize=20, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_aq[label]
        title = f"AQ: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
        ax.set_title(title, fontsize=7, color='black')

    results_snr_aq_220, = get_snrs(gw_event, start_time, gr220=False, aq220=True, gr221=False, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='red', label='AQ run with high settings and 220 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_220[0]:.2f}$_{{-{results_snr_aq_220[1]:.2f}}}^{{+{results_snr_aq_220[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)
    
    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/full_posts/GW_{gw_event}_220_AQ.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/full_posts/GW_{gw_event}_220_AQ.png', dpi=150)
        
    plt.close(figure)

def single_corner_plot_GR_221(gw_event, dir_gr_221, gw_catalogue_flag, mf_truth, start_time):

    data_gr = np.loadtxt(dir_gr_221, skiprows=1)

    data_gr_physical = data_gr[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'A2221', 'phi2220', 'phi2221']

    results_gr = get_median_and_errors(data_gr_physical, labels)

    truths = [None, mf_truth, None, None, None, None, None, None]

    figure = corner.corner(data_gr_physical, labels=labels, color='blue',
                        truths = truths,
                        truth_color='k',  contourf_kwargs={'colors': ramp('blue')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} GR 221 mode at {start_time}', fontsize=20, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_gr[label]
        title = f"GR: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
        ax.set_title(title, fontsize=7, color='black')
    
    results_snr_gr_221, = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=True, aq221=False)

    #Legends

    plot = [mlines.Line2D([], [], color='blue', label='GR run with high settings and 221 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_GR = {results_snr_gr_221[0]:.2f}$_{{-{results_snr_gr_221[1]:.2f}}}^{{+{results_snr_gr_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)
    #saving the corner plot
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/full_posts/GW_{gw_event}_full_221_GR.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/full_posts/GW_{gw_event}_full_221_GR.png', dpi=150)
            
    plt.close(figure)

def single_corner_plot_AQ_221(gw_event, dir_aq_221, gw_catalogue_flag, mf_truth, start_time):

    data_aq = np.loadtxt(dir_aq_221, skiprows=1)

    data_aq_physical = data_aq[:, :-2]

    labels = ['psi', 'Mf', 'af', 'alpha', 'cosiota', 'A2220', 'A2221', 'phi2220', 'phi2221']

    results_aq = get_median_and_errors(data_aq_physical, labels)

    truths = [None, mf_truth, None, None, None, None, None, None, None]

    figure = corner.corner(data_aq_physical, labels=labels, color='blue',
                        truths = truths,
                        truth_color='k',  contourf_kwargs={'colors': ramp('red')}, **copy.deepcopy(STYLE))

    figure.suptitle(f'GW{gw_event} AQ 221 mode at {start_time}', fontsize=20, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_aq[label]
        title = f"AQ: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
        ax.set_title(title, fontsize=7, color='black')

    results_snr_aq_221, = get_snrs(gw_event, start_time, gr220=False, aq220=False, gr221=False, aq221=True)

    #Legends

    plot = [mlines.Line2D([], [], color='red', label='AQ run with high settings and 221 mode')]

    snr = [mlines.Line2D([], [], color='black', label=f'SNR_AQ = {results_snr_aq_221[0]:.2f}$_{{-{results_snr_aq_221[1]:.2f}}}^{{+{results_snr_aq_221[2]:.2f}}}$')]
    #figure.legend(handles=[blue_line, red_line], loc='upper left', bbox_to_anchor=(0.45, 0.95), fontsize=10)
    #figure.legend(handles = snr, loc = 'upper left', bbox_to_anchor=(0.45, 0.78), fontsize = 10)

    # Legend 1: top right corner of the figure
    plot_legend = figure.legend(handles=plot, loc='upper right',
                            bbox_to_anchor=(0.98, 0.98), fontsize=20)

    # Measure where it ended up (figure coordinates)
    figure.canvas.draw()
    bb = plot_legend.get_window_extent().transformed(figure.transFigure.inverted())

    # Legend 2: right edge aligned with legend 1, directly underneath it
    snr_legend = figure.legend(handles=snr, loc='upper right',
                           bbox_to_anchor=(bb.x1, bb.y0 - 0.01), fontsize=20)
    
    #saving the corner plot 
    if gw_catalogue_flag == 'GWTC3':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC3/full_posts/GW_{gw_event}_221_AQ.png', dpi=150)
    elif gw_catalogue_flag == 'GWTC4':
        figure.savefig(f'/groups/astro/ctj304/corner_plots/{start_time}/GWTC4/full_posts/GW_{gw_event}_221_AQ.png', dpi=150)
    
    plt.close(figure)
    
def plots_for_tests(gw_event):

    dir_1 = f'/groups/astro/ctj304/catalogue_events_prod/GW230927_043729_220_GR_TEST_NORMAL/Nested_sampler/posterior.dat' #normal
    dir_2 = f'/groups/astro/ctj304/catalogue_events_prod/GW230927_043729_220_GR_DOWNLOAD/Nested_sampler/posterior.dat' #download
    dir_3 = f'/groups/astro/ctj304/catalogue_events_prod/GW230927_043729_220_GR_TEST_DCC/Nested_sampler/posterior.dat' #dcc
    
    data_1 = np.loadtxt(dir_1, skiprows=1)
    data_2 = np.loadtxt(dir_2, skiprows=1)
    data_3 = np.loadtxt(dir_3, skiprows=1)

    data_1_physical = data_1[:, :-2]
    data_2_physical = data_2[:, :-2]
    data_3_physical = data_3[:, :-2]

    labels = ['psi', 'Mf', 'af', 'cosiota', 'A2220', 'phi2220']

    results_1 = get_median_and_errors(data_1_physical, labels)
    results_2 = get_median_and_errors(data_2_physical, labels)
    results_3 = get_median_and_errors(data_3_physical, labels)


    figure = corner.corner(data_1_physical, labels=labels, color='blue',
                         truth_color='blue')
    corner.corner(data_2_physical, labels=labels, color='red', fig=figure,
               truth_color='red')
    corner.corner(data_3_physical, labels=labels, color='green', fig=figure,
               truth_color='green')

    figure.suptitle(f'GW{gw_event} normal vs download vs dcc 220 v2', fontsize=15, y=0.8, x=0.8, color='black')

    # Manually write both models' medians/errors as titles
    ndim = len(labels)
    axes = np.array(figure.axes).reshape((ndim, ndim))
    
    #prints medians and errors
    for i, label in enumerate(labels):
        ax = axes[i, i]
        med1, minus1, plus1 = results_1[label]
        med2, minus2, plus2 = results_2[label]
        med3, minus3, plus3 = results_3[label]
        title = (f"normal: {med1:.2f}$_{{-{minus1:.2f}}}^{{+{plus1:.2f}}}$\n"
                f"download: {med2:.2f}$_{{-{minus2:.2f}}}^{{+{plus2:.2f}}}$"
                f"dcc: {med3:.2f}$_{{-{minus3:.2f}}}^{{+{plus3:.2f}}}$")
        ax.set_title(title, fontsize=7, color='black')
    
    #Legend
    red_line = mlines.Line2D([], [], color='blue', label='normal')
    green_line = mlines.Line2D([], [], color='green', label='dcc')
    blue_line = mlines.Line2D([], [], color='red', label='download')
    figure.legend(handles=[blue_line, red_line, green_line], loc='upper right', fontsize=25)
    
    #saving the corner plot
    figure.savefig(f'/groups/astro/ctj304/corner_plots/misc/GW_{gw_event}_220_normal_download_dcc.png', dpi=150)
    
    plt.close(figure)

def get_snrs(gw_event, start_time, gr220: False, aq220: False, gr221: False, aq221: False):

    requested = [
        (gr220, "220", "GR"),
        (aq220, "220", "AQ"),
        (gr221, "221", "GR"),
        (aq221, "221", "AQ")
    ]

    results=[]

    for flag, mode, run in requested:
        if flag:
            dir_snr = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{gw_event}_{mode}_{run}_{start_time}/Nested_sampler/optimal_SNR_TD.dat'
            data_snr = np.loadtxt(dir_snr, skiprows=1)
            data_snr_physical = data_snr[:, 0]
            results_snr = get_median_and_errors_single(data_snr_physical)
            results.append(results_snr)
    return tuple(results)  # Return a tuple of results for all requested modes/runs

def get_truths(gw_event, gw_catalogue_flag, start_time):

    config = configparser.ConfigParser()
    config.optionxform = str

    config.clear()

    if gw_catalogue_flag == 'GWTC3':
        file=config.read(f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time}/GWTC3/config_GW_{gw_event}_220_GR_prod_download_{start_time}.ini')
        print(file)
        mf_truth = config.getfloat("Priors", "mf-time-prior")
    elif gw_catalogue_flag == 'GWTC4':
        config.read(f'/groups/astro/ctj304/pyring/pyRing/config_files/catalogue_{start_time}/GWTC4/config_GW_{gw_event}_220_GR_prod_download_{start_time}.ini') #all models use same true values
        mf_truth = config.getfloat("Priors", "mf-time-prior")

    return mf_truth #only Mf true for now


gw_catalogue_GWTC4 = [
    '230601_224134',
    '230609_064958',
    '230628_231200',
    '230811_032116',
    '230814_061920',
    '230824_033047',
    '230914_111401',
    '230922_020344',
    '230922_040658',
    '230924_124453',
    '230927_043729',
    '230927_153832',
    '230928_215827',
    '231001_140220',
    '231028_153006',
    '231102_071736',
    '231108_125142',
    '231213_111417',
    '231223_032836',
    '231226_101520',
]


gw_catalogue_GWTC3 = [
    #'150914',
    #'170104',
    #'170814',
    #'170823',
    #'190408_181802',
    #'190512_180714',
    #'190513_205428',
    #'190519_153544',
    #'190521_030229',
    #'190521_074359',
    #'190602_175927',
    #'190706_222641',
    #'190708_232457',
    #'190727_060333',
    #'190828_063405',
    '190910_112807',
    #'190915_235702',
    #'191109_010717',
    #'191222_033537',
    #'200129_065458',
    #'200224_222234',
    #'200311_115853',
]

gw_catalogue_flag = 'GWTC3'
gw_catalogue = gw_catalogue_GWTC3
prev_date, k = None, 0
start_time = '10M'
if gw_catalogue is gw_catalogue_GWTC3:
    gw_catalogue_flag = 'GWTC3'
    for i in gw_catalogue:
        dir_aq_220 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_220_AQ_{start_time}/Nested_sampler/posterior.dat'
        dir_gr_220 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_220_GR_{start_time}/Nested_sampler/posterior.dat'
        dir_aq_221 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_221_AQ_{start_time}/Nested_sampler/posterior.dat'
        dir_gr_221 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_221_GR_{start_time}/Nested_sampler/posterior.dat'
        directory_gwtc = f'/groups/astro/ctj304/GWTC3/rin/'
        
        date = i[:6]
        k = k+1 if date == prev_date else 0
        prev_date = date 

        dir_gwtc = sorted(glob.glob(f'{directory_gwtc}rin_S{date}*_pyring_Kerr_220_10M.h5'), key = suffix_key)

        if k<len(dir_gwtc): dir_gwtc = [dir_gwtc[k]]

        print(f'name: {i}, directory: {dir_gwtc}')
        #out = get_snrs(i, gr220=True, aq220=False, gr221=False, aq221=False)
        #print(type(out), len(out))
        #print(out)
        mf_truth = get_truths(i, gw_catalogue_flag = gw_catalogue_flag, start_time=start_time)
        #corner_plots_220(i, dir_aq_220,dir_gr_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #corner_plots_221(i, dir_aq_221, dir_gr_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        corner_plots_GWTC4_220(i, dir_aq_220,dir_gr_220,dir_gwtc, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #corner_plots_GWTC4_221(i, dir_aq_221,dir_gr_221,dir_gwtc, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #cp_mfvsalpha_220(i, dir_aq_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #cp_mfvsalpha_221(i, dir_aq_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #corner_plots_GRvsGWTC_220(i, dir_gr_220, dir_gwtc, gw_catalogue_flag, mf_truth=mf_truth, start_time=start_time)
        #corner_plots_GRvsGWTC_221(i, dir_gr_221, dir_gwtc, gw_catalogue_flag, mf_truth=mf_truth, start_time=start_time)
        #posterior_alpha_220(i, dir_aq_220, gw_catalogue_flag, start_time=start_time)
        #posterior_alpha_221(i, dir_aq_221, gw_catalogue_flag, start_time=start_time)
        #corner_plots_GWTC_only(i, dir_gwtc, gw_catalogue_flag, mf_truth=mf_truth, start_time=start_time) 
        #single_corner_plot_GR_220(i, dir_gr_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #single_corner_plot_AQ_220(i, dir_aq_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #single_corner_plot_GR_221(i, dir_gr_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #single_corner_plot_AQ_221(i, dir_aq_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #plots_for_tests(i)
elif gw_catalogue is gw_catalogue_GWTC4:
    gw_catalogue_flag = 'GWTC4'
    for i in gw_catalogue:
        dir_aq_220 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_220_AQ_{start_time}/Nested_sampler/posterior.dat'  
        dir_gr_220 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_220_GR_{start_time}/Nested_sampler/posterior.dat'
        dir_aq_221 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_221_AQ_{start_time}/Nested_sampler/posterior.dat'
        dir_gr_221 = f'/groups/astro/ctj304/catalogue_events_prod/catalogue_events_prod_{start_time}/GW{i}_221_GR_{start_time}/Nested_sampler/posterior.dat'   
        dir_gwtc =  f'/groups/astro/ctj304/GWTC4/pyring/samples/GW{i}/GW{i}_Kerr_220_np_10M/pyRing_GW{i}_Kerr_220_np_10M.dat'
        mf_truth = get_truths(i, gw_catalogue_flag = gw_catalogue_flag, start_time=start_time)
        #out = get_snrs(i, gr220=True, aq220=False, gr221=False, aq221=False)
        #print(type(out), len(out))
        #print(out)
        corner_plots_220(i, dir_aq_220,dir_gr_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        corner_plots_221(i, dir_aq_221, dir_gr_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        corner_plots_GWTC4_220(i, dir_aq_220,dir_gr_220,dir_gwtc, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        corner_plots_GWTC4_221(i, dir_aq_221,dir_gr_221,dir_gwtc, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        cp_mfvsalpha_220(i, dir_aq_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        cp_mfvsalpha_221(i, dir_aq_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        corner_plots_GRvsGWTC_220(i, dir_gr_220, dir_gwtc, gw_catalogue_flag, mf_truth= mf_truth, start_time=start_time)
        corner_plots_GRvsGWTC_221(i, dir_gr_221, dir_gwtc, gw_catalogue_flag, mf_truth=mf_truth, start_time=start_time)
        posterior_alpha_220(i, dir_aq_220, gw_catalogue_flag, start_time=start_time)
        posterior_alpha_221(i, dir_aq_221, gw_catalogue_flag, start_time=start_time)
        corner_plots_GWTC_only(i, dir_gwtc, gw_catalogue_flag, mf_truth=mf_truth, start_time=start_time) 
        single_corner_plot_GR_220(i, dir_gr_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        single_corner_plot_AQ_220(i, dir_aq_220, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        single_corner_plot_GR_221(i, dir_gr_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        single_corner_plot_AQ_221(i, dir_aq_221, gw_catalogue_flag, mf_truth = mf_truth, start_time=start_time)
        #plots_for_tests(i)

