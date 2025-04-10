import sys
import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from operator import sub
import datetime
import re 
import matplotlib.cm as cm  # Import colormap module

datafiles = {
    "det.txt": "Determinant Error", 
    "sjac.txt": "Scaled Jacobian Error",
    "int_err.txt": "Integrability Error",
    "balign.txt": "Boundary Alignment Error",
    "combed_smooth_grad.txt": "Combed Dirichlet of Jacobian",
    "inv_elements_fraction.txt": "Fraction of inverted elements",
    "aniso.txt": "Anisotropy Error"
    }


# Reads data from one file and calculates the mean and standard deviation.
# Adjusts determinant and scaled Jacobian values to measure error (away from 1.0)
def meshStats(datafile):
    # print(datafile)
    data = np.loadtxt(datafile)
    cleaneddata = data[~np.isnan(data)]
    count = cleaneddata.size
    
    # fix some of the data files
    fname = os.path.basename(datafile)
    if fname == "det.txt":
        a = np.abs(cleaneddata)
        det_mean = a.mean()
        a = a * (1. / det_mean)
        d = np.abs(1 - a)
#     elif fname == "inv_elements_fraction.txt":
#         d = cleaneddata
#         if d.mean() > 0.5:
#  # error 
#             print("bleep")
#             print(cleaneddata)
#             print(datafile)
#             print("bleep")
#             raise Exception("Invalid data file type: " + datafile)

 
 #           print("bleep")
 #           print(cleaneddata)
    #     print("bleep")
    #     print(cleaneddata)

    elif fname == "aniso.txt":
        d = np.abs(1 - cleaneddata)
    elif fname == "sjac.txt":
        d = np.abs(1 - cleaneddata)
    elif fname == "int_err.txt" or fname == "balign.txt" or fname == "combed_smooth_grad.txt" or fname == "inv_elements_fraction.txt":
        d = cleaneddata
    else:
        raise Exception("Invalid data file type: " + datafile)
    
    mean = d.mean()
    std = d.std()
    percent10 = np.percentile(d, 10)
    percent90 = np.percentile(d, 90)
    p10error = max(0.0, mean - percent10)
    p90error = max(0.0, percent90 - mean)
    return (mean.item(), std.item(), p10error, p90error)

# Recursively searches a folder for data files and collects per-mesh stats.
# Return a dictionary {stat-name -> {mesh-name -> {method-name -> (mean, stddev, ...)}}}
def collectMeshStats(basefolder):
    stats = defaultdict(lambda: defaultdict(dict))
    
    methods = os.listdir(folder)

    for method in methods:
        methodfolder = os.path.join(folder, method)
        for example in os.listdir(methodfolder):
            examplefolder = os.path.join(methodfolder, example)
            if os.path.isdir(examplefolder):
                found_inv_elem_count = False
                for file in os.listdir(examplefolder):
                                            # # compute number of inverted elements in parameterization
                    if file in datafiles:
                        filepath = os.path.join(examplefolder, file)
                        result = meshStats(filepath)                        
                        stats[file][example][method] = result
                            # a = np.abs(cleaneddata)
                            # det_mean = a.mean()
                            # a = a * (1. / det_mean)
                            # d = np.abs(1 - a)
                            # if fname == "det.txt":

                if not found_inv_elem_count:
                # if True:    
                    fname = "inv_elements_fraction.txt"
                    filepath = os.path.join(examplefolder, fname)

                    has_inv_elem_count = os.path.isfile(filepath)
                    if not has_inv_elem_count:
                    # if True:
                        detpath = os.path.join(examplefolder, "det.txt")
                        if os.path.isfile(detpath):
                            data = np.loadtxt(detpath)
                            cleaneddata = data[~np.isnan(data)]
                            sgn_flip = np.sum(cleaneddata < 0)
                            data_sorted = np.sort(cleaneddata)
                            if (data_sorted[len(cleaneddata) // 2] < 0):
                                sgn_flip = len(cleaneddata) - sgn_flip
                            with open(os.path.join(examplefolder, "inv_elements_fraction.txt"), 'w') as f:
                                print(sgn_flip / len(cleaneddata))
                                f.write(str(sgn_flip / len(cleaneddata)))
                            stats["inv_elements_fraction"][example][method] = (sgn_flip / len(cleaneddata), 0,0,0)
    return stats
    
# Processes a dictionary {mesh-name -> {method-name -> (mean, stddev, ...)}} to turn the list of examples into a data series that can be plotted.
# Returns a dictionary {method-name -> ((mean1, stddev1, ...), (mean2, stddev2, ...), ...)}
# Each method's value is a list of the same size (the total number of meshes).
# Throws if some meshes don't have data for all methods.
def aggregateOneStat(meshdict):
    methods = set()
    for e in meshdict:
        for m in meshdict[e]:
            methods.add(m)
    
    methlist = list(methods)
    data = {}
    for m in methlist:
        data[m] = []
        
    problems = [];
    
    for e in meshdict:
        exampleok = True
        for m in methlist:
            if m not in meshdict[e]:
                problems.append((e,m))
                exampleok = False
        if exampleok:
            for m in methlist:
                data[m].append(meshdict[e][m])
                        
    if problems:
        print("Warning: data is inconsistent!")
        for p in problems:
            print("Example " + p[0] + " is missing method " + p[1])        
    return data
    
# Given a dictionary {method-name -> ((mean1, stddev1, ...), (mean2, stddev2, ...), ...)}, sorts all data series based on the means of the given key method.
def sortOneStat(methodsdata, key):
    if key not in methodsdata:
        raise Exception(key + " not in data")
    values = []
    names = []
    names.append(key)
    values.append(methodsdata[key])
    othermethods = ()
    for m in methodsdata:
        if m != key:
            names.append(m)
            values.append(methodsdata[m])
    
    newvalues = list(zip(*sorted(zip(*values))))
    result = {}
    for i in range(len(names)):
        result[names[i]] = newvalues[i]
    return result
    
def sortOneStatVsRef(methodsdata, key, baseline):
    if key not in methodsdata:
        raise Exception(key + " not in data")
    values = []
    names = []
    names.append("metricguided_ray_vs_mint_diff")

    diff = [tuple(map(sub, methodsdata[key][i][0:4], methodsdata[baseline][i][0:4])) for i in range(len(methodsdata[key]))]

    values.append(diff)
    othermethods = ()
    for m in methodsdata:
        # if m != key:
        names.append(m)
        values.append(methodsdata[m])

        # values.append([ methodsdata[key][i][0:4] for i in range(len(methodsdata[m]))])
    
    newvalues = list(zip(*sorted(zip(*values))))
    result = {}
    for i in range(len(names)):
        result[names[i]] = newvalues[i]
    return result
    

def plotOneStat(methodsdata):
    fig = plt.figure(figsize=(15, 15))
    ax = plt.subplot(111)

    min_x, max_x = float('inf'), float('-inf')
    min_y, max_y = float('inf'), float('-inf')

    sorted_methods = sorted(methodsdata.keys())

    cmap = cm.get_cmap('tab20')  # Or any other colormap
    num_methods = len(sorted_methods)
    colors = []
    for i in range(num_methods):
        color = cmap((i+2) * 1./ 20.)
        colors.append(color)
        # ax.plot([], [], color=color, label=sorted_methods[i])
    # colors = [cmap((i + .8) - 1 ) for i in np.linspace(0, 1, num_methods)]

    for method_index, method in enumerate(sorted_methods):
        data = methodsdata[method]
        mean, stddev, p10, p90 = zip(*data)
        style = '-'
        alpha = 1.0
        linewidth=2

        if not ("mint" in method or "metricguided" in method):
            style = ':'
            alpha = .9
            linewidth=1.

        if "metricguided" in method:
            alpha = 0.90

        if "noint" in method:
            # style = ':'
            linewidth=2.5
            # alpha = 0.7

        if "diff" in method:
            style = 'dashdot'
            # alpha = 0.60
            alpha = 1
            linewidth=1.5

        if "lmff" in method:
            style = '-.'

        color = colors[method_index]
        ax.plot(range(len(mean)), mean, label=method, linestyle=style, alpha=alpha, color=color, linewidth=linewidth * 3)

        if "mint" in method or "metricguided" in method or "diff" in method:
            min_x = min(min_x, 0)
            max_x = max(max_x, len(mean) - 1)
            min_y = min(min_y, min(mean))
            max_y = max(max_y, max(mean))

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    ax.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t)) # Sort legend entries

    legend_fig, legend_ax = plt.subplots(figsize=(8, 1))  # Adjust size as needed
    legend = legend_fig.legend(handles, labels, loc='center',  # Center on the legend figure
                              ncol=4,  # Adjust columns
                              frameon=False)  # Remove frame

    legend_ax.axis('off')  # Turn off axes for the legend figure
    legend_fig.tight_layout() # Remove extra whitespace around the legend

    # Save the legend to a separate file (assuming plots_folder is defined elsewhere)
    legend_fig.savefig(os.path.join(plots_folder, "legend.png"), dpi=300, bbox_inches='tight')  
    plt.close(legend_fig) # Close the legend figure

    # ax.get_legend().remove()  # Remove the legend from the main plot
    plt.axhline(y=0, color='black', linestyle='-')

    return fig  # Return the figure object if you need it later

    # handles, labels = ax.get_legend_handles_labels()
    # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    # ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))  


# # Plots curves with error bars given data in a dictionary {method-name -> ((mean1, stddev1, ...), (mean2, stddev2, ...), ...)}
# def plotOneStat(methodsdata):
#     fig = plt.figure(figsize=(20,6))







# Extract flags from ourmethod
def extract_flags(ourmethod_str):
    flags = {}
    match_o = re.search(r"o_([\d\.eE\+\-]+)", ourmethod_str)
    match_u = re.search(r"u_([\d\.eE\+\-]+)", ourmethod_str)
    match_feat = re.search(r"feat_([\d\.eE\+\-]+)", ourmethod_str)

    if match_o:
        flags['o'] = match_o.group(1)
    if match_u:
        flags['u'] = match_u.group(1)
    if match_feat:
        flags['feat'] = match_feat.group(1)

    return flags


# ourmethod = "mint_u_1e-4_o_1e-1_feat_1e+3_seamless"
# ourmethod = "mint_o_5e-1_u_1e-2_crs_5e1+hard_seamless"
# ourmethod = "mint_o_1e2_u_2.5_hard_seamless"
# ourmethod = "mint_o_1e0_u_2.5e-2_crs_hard_seamless"
# ourmethod = "mint_o_1e-1_u1e-4_seamless"
# ourmethod = "mint_o_1e0_u_2.5e-2_crs_hard_seamless"
# ourmethod = "mint_o_5e-1_u_1e-2_seamless" # pretty good 
# ourmethod = "mint_0_5e-1_u_1e-2_feat_1e4_maybe_seamless" 
# ourmethod = "mint_o_5e-1_u_1e-2_feat_5e1_seamless"
# ourmethod = "mint_o_1e2_u_2.5_crs_1e3_seamless"
# ourmethod = "mint_o_5e-1_u_1e-2_crs_1e4_seamless"
# ourmethod = "mint_o_5e-1_u_1e-2_feat_1e+4_seamless"
# ourmethod = "mint_o_1e-1_u_1e-3_crs_1e3_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e-01_u_1.00e-02_vsc_1.00e-09_feat_1.00e+03_bsdf_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_2.50e-01_u_1.00e-03_vsc_1.00e-09_feat_1.00e+03_bsdf_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e+01_u_1.00e-03_vsc_1.00e-09_feat_1.00e+04_bsdf_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_7.50e-02_u_1.00e-03_vsc_1.00e-09_feat_1.00e+03_bsdf_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_7.50e-02_u_1.00e-03_vsc_1.00e-09_feat_1.00e+03_bsdf_seamless"
# ourmethod = "o_2.50e+01_u_1.00e-04_vsc_1.00e-09_feat_1.00e+04_seamless"

# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_2.50e-01_u_5.00e-03_vsc_1.00e-09_bsdf_seamless"
# baseline = "metricguided_solids_intgrid"
# baseline = "mint_mesh"
# baseline = "lmff_solids_ray_seamless"
# baseline = "mint_o_1e-1_u_1e-3_crs_1e3_seamless"
# baseline = "arff_solids_Octa_mMBO_seamless"
# baseline = "mint_mesh_solids_final_seamless"
# baseline = "metricguided_seamless"
# baseline = "metricguided_intgrid"

# baseline = "metricguided_ray"
# baseline = "mint_mesh_solids_noint_o_0.5_u_1e-5_p_1e-6_seamless"

# baseline = "lmff_solids_locally_meshable_intgrid"
# baseline = "mint_mesh_solids_noint_o_1e-1_u_1e-3_p_1e-6_intgrid"
baseline = "mint_mesh_o_0.5_u_1e-5_p_1e-6_noint"



# ourmethod ="mint_not_sure_what_else_o_1.00e+00_u_1.00e-04_vsc_1.00e-09_feat_1.00e+04_bsdf_seamless"; # o_1.00e-01_u_1.00e-04_vsc_1.00e-09_feat_5.00e+01
# ourmethod = "mint_o_1.00e-01_u_1.00e-04_vsc_1.00e-09_feat_5.00e+01_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e+00_u_1.00e-04_vsc_1.00e-09_feat_1.00e+04_bsdf_seamless"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e+00_u_1.00e-04_vsc_1.00e-09_feat_1.00e+04_bsdf_intgrid"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e-01_u_1.00e-03_vsc_1.00e-09_feat_1.00e+03_bsdf_intgrid"
# ourmethod = "mint_ss_1.00e+00_is_1.00e-08_o_1.00e-01_u_1.00e-03_vsc_1.00e-09_feat_1.00e-01_bsdf_seamless"
# ourmethod = "mint_mesh_with_crease_regularity_term_intgrid"
# ourmethod = "mint_o_1e-1_u_1e-3_crs_1e3_intgrid"
# ourmethod = "mint_plane_w_feature_align_seamless"
# ourmethod = "mint_plane_1e-3_seamless"
# ourmethod = "mint_o_1e-1_u_1e-3_seamless"
# ourmethod = "mint_mesh_plane_0_crease5e1_converged_seamless"
# ourmethod = "mint_mesh_plane_1e-3_crease_0_seamless"

# ourmethod = "mint_mesh_polyhedral_plane_1e-5_seamless"
# ourmethod ="mint_mesh_plane_1e-2_intgrid"
# ourmethod = "mint_mesh_plane_1e-3_crease_0_intgrid"
# ourmethod = "mint_mesh_solids_o_1e-1_u_1e-3_p_1e-6_seamless"
# ourmethod = "mint_mesh_solids_rerun_seamless"
# ourmethod = "mint_mesh_intgrid"

# ourmethod =  "mint_mesh_solids_o_7.5e-1_u_1e-6_p_1e-7_seamless"
# ourmethod = "mint_mesh_solids_o_0.5_u_1e-5_p_1e-6_intgrid"
# ourmethod = "mint_mesh_solids_o_0.5_u_1e-5_p_1e-6_seamless"
# ourmethod = "mint_mesh_solids_o_0.5_u_1e-5_p_1e-6_rerun_seamless"

ourmethod = "mint_mesh_o_0.5_u_1e-5_p_1e-6"



folder = sys.argv[1]
data = collectMeshStats(folder)

######
# save to  file 
######

flags = extract_flags(ourmethod)
flags_str = ""
if 'o' in flags:
    flags_str += f"o_{flags['o']}_"
if 'u' in flags:
    flags_str += f"u_{flags['u']}_"
if 'feat' in flags:
    flags_str += f"feat_{flags['feat']}_"
flags_str = flags_str.rstrip("_") # remove the trailing underscore
    
# Create the plots directory with timestamp and flags
now = datetime.datetime.now()
timestamp_str = now.strftime("%Y-%m-%d_%H-%M-%S")
plots_folder = os.path.join(folder + "/plots", f"{timestamp_str}_{flags_str}")
os.makedirs(plots_folder, exist_ok=True)


######
# run the plotting
######

    

for d in data:
    agd = aggregateOneStat(data[d])
    # sortd = sortOneStat(agd,ourmethod)
    sortd = sortOneStatVsRef(agd,ourmethod,baseline)
    plotOneStat(sortd)
    plt.ylabel(datafiles[d])

    # Save the figure in the timestamped plots folder
    fig_filename = os.path.splitext(d)[0] + ".png"
    fig_path = os.path.join(plots_folder, fig_filename)
    plt.savefig(fig_path, dpi=300)
    plt.clf()