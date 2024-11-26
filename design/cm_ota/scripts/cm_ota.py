
import sys, os, getpass, shutil, operator, collections, copy, re, math
import matplotlib.ticker as mticker
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm
ROAR_HOME = os.environ["ROAR_HOME"]
ROAR_LIB = os.environ["ROAR_LIB"]
ROAR_SRC = os.environ["ROAR_SRC"]
ROAR_CHARACTERIZATION = os.environ["ROAR_CHARACTERIZATION"]
ROAR_DESIGN = os.environ["ROAR_DESIGN"]
sys.path.append(ROAR_SRC)
from cid import *
from matplotlib.font_manager import FontProperties
plt.rcParams['svg.fonttype'] = 'none'

def parallel(x1, x2):
    return 1/((1/x1) + 1/x2)

def create_spice_netlist_for_place_and_route(nf1_2, nf3_4, nf5_6, nf7_8, nf9_10):
    print("TODO")

def total_current_ota_v2(ncorner, pcorner, alpha, gbw, cload, kgm_n, kgm_p, gain_spec, thermal_noise_spec):
    f1 = 2*math.pi*gbw
    p_factor = 1
    k = 1.380649e-23
    T = 300.15
    gamma = 2/3
    gamma = 8/3
    kgm1 = kgm_n
    kgm2 = kgm_n
    kgm5 = kgm_n
    kgm6 = kgm_n
    kgm3 = kgm_p
    kgm4 = kgm_p
    kgm7 = kgm_p
    kgm8 = kgm_p
    kcgs_1 = ncorner.lookup(param1="kgm", param2="kcgs", param1_val = kgm1)
    kcds_1 = ncorner.lookup(param1="kgm", param2="kcds", param1_val = kgm1)
    kcgs_8 = pcorner.lookup(param1="kgm", param2="kcgs", param1_val=kgm8)
    kcgd_8 = pcorner.lookup(param1="kgm", param2="kcgd", param1_val=kgm8)
    kcds_8 = pcorner.lookup(param1="kgm", param2="kcds", param1_val=kgm8)
    kcgs_6 = ncorner.lookup(param1="kgm", param2="kcgs", param1_val=kgm6)
    kcgd_6 = ncorner.lookup(param1="kgm", param2="kcgd", param1_val=kgm6)
    kcds_6 = ncorner.lookup(param1="kgm", param2="kcds", param1_val=kgm6)
    kcgs_4 = pcorner.lookup(param1="kgm", param2="kcgs", param1_val=kgm4)
    kcds_4 = pcorner.lookup(param1="kgm", param2="kcds", param1_val=kgm4)
    kgds_6 = ncorner.lookup(param1="kgm", param2="kgds", param1_val=kgm6)
    kgds_8 = ncorner.lookup(param1="kgm", param2="kgds", param1_val=kgm8)
    gain = kgm1/(kgds_6 + kgds_8)
    kcout = kcgd_8 + kcds_8 + kcgd_6 + kcds_6
    beta_num = kgm4 - 2*math.pi*alpha*gbw*kcgs_4
    beta_denom = (2*math.pi*alpha*gbw*(kcgs_8 + kcgd_8))
    beta = beta_num/beta_denom
    #m1_current_term1 = 2*math.pi*gbw*cload/(beta*kgm1)
    #m1_current_term2 = 1/(1 - 2*math.pi*gbw*kcout/kgm1)
    #m1_current = m1_current_term1*m1_current_term2
    #total_current = m1_current + beta*m1_current
    kcs_8 = kcgs_8 + kcgd_8
    kcs_6 = kcds_6 + kcgs_6
    total_current_num = alpha*cload*f1*f1*(kcs_8 - kcgs_4) + f1*cload*kgm4
    total_current_denom = (kgm4 - alpha*f1*kcgs_4)*(kgm1 - f1*kcout)
    total_current = total_current_num/total_current_denom
    m1_current = total_current/(1 + beta)
    m8_current = total_current - m1_current
    m6_current = m8_current
    m4_current = m1_current
    m2_cload = (kcgs_1 + kcds_1)*m1_current + kcs_8*m8_current + kcds_4*m1_current
    m4_cload = m2_cload
    m6_cload = ((kcs_6 + kcds_8 + kcgs_8)*m6_current) + cload
    m8_cload = m6_cload
    m1_gm = kgm1*m1_current
    m6_gm = kgm2*m6_current
    m8_gm = m6_gm
    m4_gm = m1_gm
    ft_m8 = m8_gm/(2*math.pi*m8_cload)
    ft_m2 = m1_gm/(2*math.pi*m2_cload)
    ft_m4 = m4_gm/(2*math.pi*m4_cload)
    ft_m6 = m6_gm/(2*math.pi*m6_cload)

    thermal_rms_noise = k*T*(ft_m2/(kgm1*m1_current) + ft_m4/(kgm1*m1_current) + ft_m6/(kgm2*m6_current) + ft_m8/(kgm2*m6_current))
    #thermal_rms_noise = (gamma*k*T*ft_m1)/(m1_current*kgm1)
    #thermal_rms_noise = 2*((8/3)*k*T)*(1/m1_gm)*(1 + 1 + 1/beta + m6_gm/(beta*beta + m1_gm))
    beta_valid = True
    gain_valid = True
    thermal_noise_valid = True
    if beta < 1:
        beta_valid = False
    if gain < gain_spec:
        gain_valid = False
    if thermal_rms_noise < thermal_noise_spec:
        thermal_noise_valid = False
    return total_current, m1_current, m6_current, beta, kcout, gain, thermal_rms_noise, beta_valid, gain_valid, thermal_noise_valid, kcout


def plot_results_krummenechar_ota_stage1(nom_ncorner, nom_pcorner, alpha, gain, bw, cload, therm_noise, fig, ax1, ax2, ax3, ax4,
                                         color_map, beta_color, gain_color, alpha_graph, alpha_region,hatch_mark, marker_size, line_style,
                                         kgm_n_max, kgm_p_max, map_label, edge_color, legend_str):
    gbw = gain*bw
    kgm_n_v = nom_ncorner.df["kgm"]
    kgm_p_v = nom_pcorner.df["kgm"]
    max_n = max(kgm_n_v)
    max_p = max(kgm_p_v)
    min_n = min(kgm_n_v)
    min_p = min(kgm_p_v)
    kgm_min = min(min_n, min_p)
    kgm_max = max(max_n, max_p)
    current_max = 5000
    current_min = 120
    kgm_min = 0.1
    num_samples = 20
    if kgm_min < 0:
        kgm_min = 0.001
    kgm_min = 0.1
    kgm_max = 26.5
    kgm_vals = np.linspace(kgm_min, kgm_max, num_samples)
    kgm1_grid, kgm2_grid = np.meshgrid(kgm_vals, kgm_vals)
    z = np.zeros_like(kgm1_grid)
    z_log = np.zeros_like(kgm1_grid)
    kcout = np.zeros_like(kgm1_grid)
    kcout_log = np.zeros_like(kgm1_grid)
    gain_v_v = np.zeros_like(kgm1_grid)
    gain_db = np.zeros_like(kgm1_grid)
    beta = np.zeros_like(kgm1_grid)
    beta_valid_grid = np.zeros_like(kgm1_grid)
    gain_valid_grid = np.zeros_like(kgm1_grid)
    kco_max = 1e-05
    if map_label:
        kgm_init_n_max = kgm_n_max
        kgm_init_p_max = kgm_p_max

    for i in range(len(kgm_vals)):
        for j in range(len(kgm_vals)):
            total_current, m1_current_i_j, m6_current_i_j, beta_i_j, kcout_i_j, gain_i_j, thermal_rms_noise_i_j, beta_valid_i_j, gain_valid, thermal_noise_valid, kc_out= total_current_ota_v2(nom_ncorner, nom_pcorner, alpha, gbw, cload,
                                                                                                                           kgm_vals[i], kgm_vals[j], gain_spec=gain,
                                                                                                                           thermal_noise_spec=therm_noise)
            if kgm_vals[j] > kgm_p_max:
                total_current = -1
                #kcout_i_j = np.nan
            if kgm_vals[i] > kgm_n_max:
                total_current = -1
                #kcout_i_j = np.nan
            total_current = total_current * 1e6
            gain_db_i_j = 20*math.log10(gain_i_j)
            if total_current < current_min:
                total_current = np.nan
                gain_i_j = np.nan
                gain_db_i_j = np.nan
                gain_valid = False
            if total_current > current_max:
                total_current = current_max
            if kcout_i_j > kco_max or kgm_vals[j] > 18.8 or kgm_vals[i] > 26.26:
                kcout_i_j = np.nan

            z[i, j] = total_current
            z_log[i, j] = math.log10(total_current)
            kcout[i, j] = kcout_i_j
            kcout_log[i, j] = math.log10(kcout_i_j)

            gain_valid_grid[i, j] = gain_valid
            gain_v_v[i, j] = gain_i_j
            gain_db[i, j] = gain_db_i_j
            beta_valid_grid[i, j] = beta_valid_i_j

    arial_bold = FontProperties(fname=ROAR_HOME + "/fonts/ArialNarrow/arialnarrow_bold.ttf")
    font_size = 12
    current_vals = [kgm_vals, z]
    gain_vals = [kgm_vals, gain_v_v]
    kco_vals = [kgm_vals, kcout]
    #ORIGINAL SCRIPT
    current_ticks = np.array([100, 250, 500, 1000, 2500, 5000])
    surf_i_total = ax1.plot_surface(kgm1_grid, kgm2_grid, z_log, lw=0.5, cmap=color_map, edgecolor=edge_color, rstride=3, cstride=3, alpha=alpha_graph, label=legend_str)
    if map_label:
        ax1.set_xlabel(r'$\mathrm{\mathcal{G}_{P}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
        ax1.set_ylabel(r'$\mathrm{\mathcal{G}_{N}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
        ax1.set_zlabel("Total Current [uA]", font=arial_bold, fontsize=font_size)
        ax1.xaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
        ax1.yaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
        ax1.zaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
        ax1.set_zticks(np.log10(current_ticks))
        ax1.set_zticklabels(current_ticks)
        ax1.set_xlim(kgm_min, kgm_p_max)
        ax1.set_ylim(kgm_min, kgm_n_max)
        ax1.set_zlim(2, 3.9)

    surf_kco_total = ax3.plot_surface(kgm1_grid, kgm2_grid, kcout_log, lw=0.5, cmap=color_map, edgecolor=edge_color, rstride=3, cstride=3, alpha=alpha_graph)
    def log_tick_formatter(val, pos=None):
        exponent = int(np.log10(val))
        return f"$10^{{{exponent}}}$" # remove int() if you don't use MaxNLocator
        # return f"{10**val:.2e}"      # e-Notation

    ax3.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    ax3.set_xlabel(r'$\mathrm{\mathcal{G}_{P}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    ax3.set_ylabel(r'$\mathrm{\mathcal{G}_{N}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    ax3.set_zlabel(r'$\mathrm{\mathcal{C}_{O}}}$ [$\mathrm{F/A}$]', font=arial_bold, fontsize=font_size)

    kco_ticks = np.array([1e-11, 1e-9, 1e-7, 1e-05])
    ax3.set_zticks(np.log10(kco_ticks))
    ax3.set_zticklabels(kco_ticks)
    if map_label:
        ax3.set_xlim(kgm_min, kgm_p_max)
        ax3.set_ylim(kgm_min, kgm_n_max)
    ax3.xaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    ax3.yaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    ax3.zaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    surf_gain_total = ax4.plot_surface(kgm1_grid, kgm2_grid, gain_v_v, lw=0.5, cmap=color_map, edgecolor=edge_color, rstride=3, cstride=3, alpha=alpha_graph)
    def log_tick_formatter(val, pos=None):
        return f"$10^{{{int(val)}}}$"  # remove int() if you don't use MaxNLocator
        # return f"{10**val:.2e}"      # e-Notation
    ax4.set_xlabel(r'$\mathrm{\mathcal{G}_{P}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    ax4.set_ylabel(r'$\mathrm{\mathcal{G}_{N}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    ax4.set_zlabel(r'$\mathrm{A_{V}}}$ [$\mathrm{V/V}$]', font=arial_bold, fontsize=font_size)
    ax4.xaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    ax4.yaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    ax4.zaxis._axinfo['grid'].update(color='gray', linestyle='--', linewidth=0.5)
    if map_label:
        ax4.set_xlim(kgm_min, kgm_p_max)
        ax4.set_ylim(kgm_min, kgm_n_max)
    # Overlay red X's for NaN values
    beta_false_mask = beta_valid_grid == False
    gain_false_mask = gain_valid_grid == False
    gain_false_mask = gain_valid_grid == False
    contour2 = ax2.contourf(kgm1_grid, kgm2_grid, z_log, levels=10, alpha=alpha_graph, cmap=color_map)
    ax2.set_xlabel(r'$\mathrm{\mathcal{G}_{P}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    ax2.set_ylabel(r'$\mathrm{\mathcal{G}_{N}}}$ [$\mathrm{V^{-1}}$]', font=arial_bold, fontsize=font_size)
    # Add gridlines for better readability
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
    cbar4 = fig.colorbar(contour2, ax=ax2, shrink=0.33)
    cbar4.set_ticks(np.log10(current_ticks))
    cbar4.set_ticklabels(current_ticks)
    if map_label:
        ax2.set_ylim(kgm_min, kgm_n_max)
        ax2.set_xlim(kgm_min, kgm_p_max)
        #ax2.invert_xaxis()
    ax2.set_ylim(0.1, 26.26)
    ax2.set_xlim(0.1, 18.8)
    cbar4.set_label(legend_str, font=arial_bold, fontsize=font_size)
    # Show the legend for the contour plot
    ax2.legend(prop=arial_bold)
    return current_vals, gain_vals, kco_vals, gain_false_mask

def cm_ota_plotting():
    arial_font = "/home/adair/Documents/CAD/roar/fonts/ArialNarrow/arialnarrow_bold.ttf"
    font_size = 6
    arial_bold = FontProperties(fname=arial_font, size=font_size)
    av= 50
    bw = 2e6
    gbw = bw * av
    therm_noise = 500e-9
    cload = 4e-12
    tan_thirty = math.tan(30*math.pi/180)
    alpha = 1/tan_thirty
    two_pi_alpha_gbw = 2*math.pi*alpha*gbw
    f2 = alpha*gbw
    nfet_nominal = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt25.csv",
                             vdd=1.8)

    pfet_nominal = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt25.csv",
                             vdd=1.8)
    nfet_hot = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt75.csv",
                             vdd=1.8)

    pfet_hot = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt75.csv",
                             vdd=1.8)
    nfet_cold = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt-25.csv",
                             vdd=1.8)

    pfet_cold = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt-25.csv",
                             vdd=1.8)

    nfet_device = CIDDevice(device_name="nfet_500n", vdd=1.8,
                            lut_directory=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500",
                            corner_list=None)
    pfet_device = CIDDevice(device_name="pfet_500n", vdd=1.8,
                            lut_directory=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500",
                            corner_list=None)
    n_list = [nfet_cold, nfet_nominal, nfet_hot]
    p_list = [pfet_cold, pfet_nominal, pfet_hot]
    fig = plt.figure(figsize=(14, 10), tight_layout=True)
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223, projection='3d')
    ax4 = fig.add_subplot(224, projection='3d')
    color_map_list = ["Blues", "Greens", "Reds"]
    beta_color_list = ["b", "g", "r"]
    line_style_list = ["r-", "g-", "b-.", "m-", "c-.", "r--", "b--", "g--", "m--"]
    alpha_graph = 0.77
    marker_size = 40
    current_vectors = []
    kgm_vectors = []
    kgm_n = [26.26, 26.26, 26.26]
    kgm_p = [18.8, 16.5, 15.25]
    edge_colors = ['royalblue', 'green', 'red']
    legend_strs = ["-25°C", "25°C", "75°C"]
    map_label = True
    current_arrays = []
    gain_arrays = []
    kco_arrays = []
    kgm_grids = []
    gain_masks = []
    alpha_region = 0.3
    hatches = ['/', '\ ', '-']
    #for i in range(len(nfet_device.corners)):
    for i in range(len(n_list)):
        nfet_corner = n_list[i]
        pfet_corner = p_list[i]
        color_map = color_map_list[i]
        beta_color = beta_color_list[0]
        gain_color = beta_color_list[0]
        edge_color = edge_colors[i]
        line_style = line_style_list[i]
        legend_str = legend_strs[i]
        hatch_mark = hatches[i]
        current, gain, kco, gain_mask = plot_results_krummenechar_ota_stage1(nom_ncorner=nfet_corner, nom_pcorner=pfet_corner, alpha=alpha,gain=av,
                                                              bw=bw, cload=cload, therm_noise=therm_noise,
                                                              fig=fig, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, color_map=color_map, beta_color=beta_color,
                                                              gain_color=gain_color, alpha_graph=alpha_graph, alpha_region=alpha_region, hatch_mark=hatch_mark,
                                                              marker_size=marker_size, line_style=line_style, kgm_n_max=kgm_n[i], kgm_p_max=kgm_p[i],
                                                              map_label=map_label, edge_color=edge_color, legend_str=legend_str)

        currents = 0
        kgms = 0
        current_vectors.append(currents)
        kgm_vectors.append(kgms)
        map_label = False
        kco_arrays.append(kco[1])
        gain_arrays.append(gain[1])
        current_arrays.append(current[1])
        kgm_grids.append(current[0])
        gain_masks.append(gain_mask)
        marker_size = marker_size - 10
        alpha_region = alpha_region - 0.1

    kgm_grid = kgm_grids[0]
    kgm_grid = np.meshgrid(kgm_grid, kgm_grid)
    convergence_grid = np.zeros((len(kgm_grid[0]), len(kgm_grid[0])))
    current0 = current_arrays[0]
    current1 = current_arrays[1]
    current2 = current_arrays[2]
    for i in range(len(kgm_grid[0])):
        for j in range(len(kgm_grid[0])):
            i0 = current0[i, j]
            i1 = current1[i, j]
            i2 = current2[i, j]
            convergence_grid[i, j] = number_convergence(i0, i1, i2, diff=0.10)
    conv_false_mask = convergence_grid == False
    ax2.contourf(kgm_grid[0], kgm_grid[1], gain_masks[0], levels=[0.5, 1], hatches=['/'], alpha=0)
    ax2.contourf(kgm_grid[0], kgm_grid[1], gain_masks[1], levels=[0.5, 1], hatches=['\ '], alpha=0)
    ax2.contourf(kgm_grid[0], kgm_grid[1], gain_masks[2], levels=[0.5, 1], hatches=['-'], alpha=0)
    ax2.contour(kgm_grid[0], kgm_grid[1], conv_false_mask, levels=[0.5, 1], color='k', alpha=1)
    ax2.set_ylim(0.1, 26.26)
    ax2.set_xlim(0.1, 18.8)
    ax2.legend(prop=arial_bold)
    return 0


def number_convergence(num1, num2, num3, diff=0.05):
    #mean = (num1 + num2 + num3)/3
    if num1 == np.nan or num2 == np.nan or num3 == np.nan:
        return False
    def within_diff(a, b):
        ratio = abs(a - b) / max(abs(a), abs(b))
        converge = ratio <= diff
        return converge
    converged = False
    converged1 = within_diff(num1, num2)
    converged2 = within_diff(num2, num3)
    converged3 = within_diff(num1, num3)
    if converged1 and converged2 and converged3:
        converged = True
    return converged


# Function to read the AC simulation data
def read_ac_simulation_data(filename):
    # Read the data file
    data = np.loadtxt(filename, skiprows=1)  # Skip the header row

    # The file format is expected to be: frequency vdb(2) vp(2)
    frequency = data[:, 0]
    magnitude_db = data[:, 4]
    phase_deg = (180/math.pi)*data[:, 6]

    return frequency, magnitude_db, phase_deg

# Function to find relevant points (DC Gain, Pole 1, Pole 2, Unity Gain Frequency)
def find_relevant_points(frequency, magnitude_db, phase_deg):
    # DC Gain: the magnitude at the lowest frequency (0 Hz or close to it)
    dc_gain = magnitude_db[0]
    # Unity Gain Frequency: the frequency where magnitude crosses 0 dB
    unity_gain_freq = None
    for i in range(len(magnitude_db) - 1):
        if magnitude_db[i] > 0 and magnitude_db[i + 1] < 0:
            unity_gain_freq = frequency[i]
            break
    # Poles: Find where phase crosses -45 degrees (Pole 1) and -135 degrees (Pole 2)
    pole_1 = None
    pole_2 = None
    for i in range(1, len(phase_deg)):
        if phase_deg[i] <= -45 and pole_1 is None:
            pole_1 = frequency[i]  # Approximate the frequency where phase crosses -45 degrees
        elif phase_deg[i] <= -135 and pole_2 is None:
            pole_2 = frequency[i]  # Approximate the frequency where phase crosses -135 degrees
            break
    return dc_gain, pole_1, pole_2, unity_gain_freq


# Function to calculate percentage difference between two values
def percentage_difference_two_values(value1, value2):
    return 100 * (value2 - value1) / value1

def percentage_difference(values):
    if not values:
        return 0
    # Calculate the mean of the values
    mean_value = sum(values) / len(values)
    # Calculate the percentage difference for each value
    percentage_differences = [(value - mean_value) / mean_value * 100 for value in values]
    # Return the average of the percentage differences
    return sum(percentage_differences) / len(percentage_differences)

# Modified function to plot magnitude and phase for ideal and extracted data
def plot_ac_results(frequency, magnitude_db, phase_deg, line_style, line_color, label_suffix="", fig=None, ax1=None, ax2=None):
    # Define font path and size within the function
    arial_font = ROAR_HOME + "/fonts/ArialNarrow/arialnarrow_bold.ttf"
    font_size = 8  # Adjust font size for labels, ticks, etc.
    font_properties = FontProperties(fname=arial_font, size=font_size)
    ideal = False
    # Create the figure and axes if not provided (for the first call)
    if fig is None or ax1 is None or ax2 is None:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.5, 2.8 * 2), dpi=300)
        #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 11), dpi=350)
        ideal = True

    # Find the relevant points
    dc_gain, pole_1, pole_2, unity_gain_freq = find_relevant_points(frequency, magnitude_db, phase_deg)
    ax1.plot(frequency, magnitude_db, linestyle=line_style, color=line_color, label=f'{label_suffix}', linewidth=1)
    # Plot the points on the graph with larger black markers
    #if dc_gain is not None:
    #    ax1.plot(100, np.interp(100, frequency, magnitude_db), 'ko', markersize=6)  # DC Gain point at 100 Hz
    #if unity_gain_freq is not None:
    #    ax1.plot(unity_gain_freq, 0, 'ko', markersize=6)  # Unity Gain point

    # Annotate the DC Gain at 100 Hz with a bubble
    magnitude_100Hz = np.interp(100, frequency, magnitude_db)
    #ax1.annotate(f'{label_suffix} DC Gain: {magnitude_100Hz:.2f} dB', xy=(100, magnitude_100Hz),
    #             xytext=(200, magnitude_100Hz + 5),
    #             fontsize=font_size, fontproperties=font_properties, ha='left',
    #             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightyellow"))

    # Annotate Unity Gain Frequency with a bubble
    #if unity_gain_freq is not None:
    #    ax1.annotate(f'{label_suffix} Unity Gain: {unity_gain_freq:.2e} Hz', xy=(unity_gain_freq, 0),
    #                 xytext=(unity_gain_freq * 1.2, -10),
    #                 fontsize=font_size, fontproperties=font_properties, ha='left',
    #                 bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightyellow"))

    # Plot phase in degrees (Ideal or Extracted)

    ax2.plot(frequency, phase_deg, linestyle=line_style, color=line_color, label=f'{label_suffix}', linewidth=1)
    # Annotate the phase margin at the unity gain frequency
    if unity_gain_freq is not None:
        phase_at_unity = np.interp(unity_gain_freq, frequency, phase_deg)  # Phase at unity gain frequency
        phase_margin = 180 + phase_at_unity  # Phase margin calculation
        #ax2.plot(unity_gain_freq, phase_at_unity, 'ko', markersize=6)  # Mark point at unity gain frequency
        #ax2.annotate(f'{label_suffix} Phase: {phase_at_unity:.1f}°\nPhase Margin: {phase_margin:.1f}°',
        #             xy=(unity_gain_freq, phase_at_unity),
        #             xytext=(unity_gain_freq * 1.2, phase_at_unity - 10),
        #             fontsize=font_size, fontproperties=font_properties, ha='left',
        #             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightyellow"))

    return fig, ax1, ax2, dc_gain, unity_gain_freq, phase_margin

# Main function to execute the script and handle comparison
def plot_spice_results():
    arial_font = ROAR_HOME + "/fonts/ArialNarrow/arialnarrow_bold.ttf"
    font_size = 8
    arial_bold = FontProperties(fname=arial_font, size=font_size)
    font_properties = arial_bold
    # Read ideal simulation data

    dc_gains = []
    unity_gains = []
    phase_margins = []

        # Read extracted simulation data
    file_neg25_sch = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/neg_25/ac_output.txt'
    freq_neg25_sch, mag_neg25_sch, phase_neg25_sch = read_ac_simulation_data(file_neg25_sch)
    # Plot extracted results and pass in the figure and axes from the ideal plot
    line_style = "-"
    line_color="blue"
    fig, ax1, ax2, dc_gain_ext, unity_gain_ext, phase_margin_ext = plot_ac_results(
        freq_neg25_sch, mag_neg25_sch, phase_neg25_sch, line_style, line_color, label_suffix="Schematic -25° C"
    )
    dc_gains.append(dc_gain_ext)
    unity_gains.append(unity_gain_ext)
    phase_margins.append(phase_margin_ext)
    file_neg25_ext = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/neg_25/ac_output_ext.txt'
    freq_neg25_ext, mag_neg25_ext, phase_neg25_ext = read_ac_simulation_data(file_neg25_ext)
    # Plot extracted results and pass in the figure and axes from the ideal plot
    line_style = "--"
    line_color="blue"
    fig, ax1, ax2, dc_gain_ext, unity_gain_ext, phase_margin_ext = plot_ac_results(
        freq_neg25_ext, mag_neg25_ext, phase_neg25_ext, line_style, line_color, label_suffix="Post Layout -25° C", fig=fig, ax1=ax1, ax2=ax2
    )
    dc_gains.append(dc_gain_ext)
    unity_gains.append(unity_gain_ext)
    phase_margins.append(phase_margin_ext)


    filename = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/25/ac_output.txt'
    frequency, magnitude_db, phase_deg = read_ac_simulation_data(filename)
    line_style = '-'
    line_color = "green"
    # Plot ideal results
    fig, ax1, ax2, dc_gain_ideal, unity_gain_ideal, phase_margin_ideal = plot_ac_results(
        frequency, magnitude_db, phase_deg, line_style, line_color, label_suffix="Schematic 25° C",  fig=fig, ax1=ax1, ax2=ax2
    )
    dc_gains.append(dc_gain_ideal)
    unity_gains.append(unity_gain_ideal)
    phase_margins.append(phase_margin_ideal)

    file_25_ext = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/25/ac_output_ext.txt'
    freq_25_ext, mag_25_ext, phase_25_ext = read_ac_simulation_data(file_25_ext)


    # Plot extracted results and pass in the figure and axes from the ideal plot
    line_style = "--"
    line_color="green"
    fig, ax1, ax2, dc_gain_ext, unity_gain_ext, phase_margin_ext = plot_ac_results(
                                 freq_25_ext, mag_25_ext, phase_25_ext, line_style, line_color, label_suffix="Post Layout 25° C", fig=fig, ax1=ax1, ax2=ax2
    )
    dc_gains.append(dc_gain_ext)
    unity_gains.append(unity_gain_ext)
    phase_margins.append(phase_margin_ext)

    file_75_ext = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/75/ac_output.txt'
    freq_75_ext, mag_75_ext, phase_75_ext = read_ac_simulation_data(file_75_ext)
    # Plot extracted results and pass in the figure and axes from the ideal plot
    line_style = "-"
    line_color="red"
    fig, ax1, ax2, dc_gain_ext, unity_gain_ext, phase_margin_ext = plot_ac_results(
        freq_75_ext, mag_75_ext, phase_75_ext, line_style, line_color, label_suffix="Schematic 75° C", fig=fig, ax1=ax1, ax2=ax2
    )
    dc_gains.append(dc_gain_ext)
    unity_gains.append(unity_gain_ext)
    phase_margins.append(phase_margin_ext)

    # Calculate percentage differences and annotate them
    dc_gain_diff = percentage_difference(dc_gains)
    unity_gain_diff = percentage_difference(unity_gains)
    phase_margin_diff = percentage_difference(phase_margins)

    file_75_ext = ROAR_DESIGN + '/cm_ota/simulation/golden_sims/75/ac_output_ext.txt'
    freq_75_ext, mag_75_ext, phase_75_ext = read_ac_simulation_data(file_75_ext)
    # Plot extracted results and pass in the figure and axes from the ideal plot
    line_style = "--"
    line_color="red"
    fig, ax1, ax2, dc_gain_ext, unity_gain_ext, phase_margin_ext = plot_ac_results(
        freq_75_ext, mag_75_ext, phase_75_ext, line_style, line_color, label_suffix="Post Layout 75° C", fig=fig, ax1=ax1, ax2=ax2
    )
    dc_gains.append(dc_gain_ext)
    unity_gains.append(unity_gain_ext)
    phase_margins.append(phase_margin_ext)
    # Calculate percentage differences and annotate them
    dc_gain_diff = percentage_difference(dc_gains)
    unity_gain_diff = percentage_difference(unity_gains)
    phase_margin_diff = percentage_difference(phase_margins)
    # Annotate percentage differences on the plots

    #ax1.annotate(f'DC Gain Diff: {dc_gain_diff:.2f}%', xy=(200, 35), fontsize=font_size, ha='left',
    #             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightcyan"))

    #ax1.annotate(f'Unity Gain Diff: {unity_gain_diff:.2f}%', xy=(1e7, 35), fontsize=font_size, ha='left',
    #             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightcyan"))

    #ax2.annotate(f'Phase Margin Diff: {phase_margin_diff:.2f}%', xy=(1e7, -90), fontsize=font_size, ha='left',
    #             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightcyan"))



    for label in ax1.get_xticklabels():
        label.set_fontproperties(arial_bold)
    for label in ax1.get_yticklabels():
        label.set_fontproperties(arial_bold)
    minor_size = font_size - 2
    ax1.set_xscale('log')
    ax1.set_ylabel('Magnitude [dB]', fontproperties=font_properties)
    ax1.grid(which='both', linestyle='--', linewidth=0.6)
    # Set the major ticks (decades)
    ax1.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
    # Add minor gridlines
    #ax1.grid(which='minor', linestyle=':', linewidth=0.4)
    # Set the minor ticks (subdivisions of decades)
    ax1.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    # Set limits for the x-axis
    ax1.set_xlim(10, 1e11)

    # Adjust tick label sizes
    ax1.tick_params(labelsize=font_size)
    #ax1.tick_params(axis='both', which='minor', labelsize=minor_size)

    ax2.set_xscale('log')
    # Set y-axis label
    ax2.set_ylabel('Phase [°]', fontproperties=font_properties)
    # Add major gridlines
    ax2.grid(which='both', linestyle='--', linewidth=0.6)
    # Set the major ticks (decades)
    ax2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
    # Add minor gridlines
    #ax2.grid(which='minor', linestyle=':', linewidth=0.4)
    # Set the minor ticks (subdivisions of decades)
    ax2.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    # Set limits for the x-axis
    ax2.set_xlim(10, 1e11)
    # Adjust tick label sizes
    ax2.tick_params(labelsize=font_size)
    #ax2.tick_params(axis='both', which='minor', labelsize=minor_size)


    #ax1.tick_params(axis='y', labelsize=12)
    for label in ax2.get_xticklabels():
        label.set_fontproperties(arial_bold)
    for label in ax2.get_yticklabels():
        label.set_fontproperties(arial_bold)
    minor_size = font_size - 2
    ax2.tick_params(axis='both', which='major', labelsize=font_size)



    #ax1.tick_params(axis='both', which='minor', labelsize=minor_size)
    ax1.legend(prop=arial_bold)
    ax2.legend(prop=arial_bold)
    # Apply tight layout and save the plots
    plt.grid(True)
    #plt.tight_layout()
    plt.savefig("ac_simulation_plot.svg", format="svg")
    plt.savefig("ac_simulation_plot.png", format="png", dpi=300)
    plt.show()


def main():

    lut_dir = ROAR_CHARACTERIZATION

    nfet_device = CIDDevice(device_name="nfet_150n", vdd=1.8,
                            lut_directory=lut_dir + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500",
                            corner_list=None)
    pfet_device = CIDDevice(device_name="pfet_150n", vdd=1.8,
                            lut_directory=lut_dir + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500",
                            corner_list=None)

    nfet_nominal = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt25.csv",
                             vdd=1.8)

    pfet_nominal = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt25.csv",
                             vdd=1.8)

    nfet_cold = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt-25.csv",
                             vdd=1.8)

    pfet_cold = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt-25.csv",
                             vdd=1.8)

    nfet_hot = CIDCorner(corner_name="nfet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/n_01v8/LUT_N_500/nfettt75.csv",
                             vdd=1.8)

    pfet_hot: object = CIDCorner(corner_name="pet_150n_nominal",
                             lut_csv=ROAR_CHARACTERIZATION + "/sky130/LUTs_SKY130/p_01v8/LUT_P_500/pfettt75.csv",
                             vdd=1.8)


    #av1 = math.sqrt(av)
    av1 = 100
    bw = 500e3
    cload1 = 50e-15
    cm_ota_plotting()
    #THIS IS GOLDEN STANDARD VALUES
    w1_2, w3_4, w5_6, w7_8 = size_ota_devices_from_kgm_and_currents(nfet_nominal, pfet_nominal,
                                                                    kgm_n=15.95843, kgm_p=5.255)
    # Test Results for smaller Kgm
    #w1_2, w3_4, w5_6, w7_8 = size_ota_devices_from_kgm_and_currents(nfet_hot, pfet_hot,
    #                                                                kgm_n=15, kgm_p=5.255)
    f1_2, f3_4, f4_5, f5_6 = get_fingers_for_align(w1_2, w3_4, w5_6, w7_8)

    plot_spice_results()
    print("Plotting Finished")
    #w1, gm1, kgm1, w2, gm2, kgm2 = krummenechar_ota_stage1(av=av1, bw=bw, cload=cload1, nfet_device=nfet_device,
    #                                                pfet_device=pfet_device, nom_ncorner=nfet_nominal, nom_pcorner=pfet_nominal)

def get_fingers_for_align(w1_2, w3_4, w5_6, w7_8):
    # Define the pitch
    pitch = 420e-9

    # Function to calculate the number of fingers
    def calculate_fingers(width):
        # Calculate the number of fingers
        num_fingers = round(width / pitch)
        # Ensure the number of fingers is even
        if num_fingers % 2 != 0:
            num_fingers += 1
        return num_fingers

    # Calculate the number of fingers for each width
    w1_2_fingers = calculate_fingers(w1_2)
    w3_4_fingers = calculate_fingers(w3_4)
    w5_6_fingers = calculate_fingers(w5_6)
    w7_8_fingers = calculate_fingers(w7_8)

    return w1_2_fingers, w3_4_fingers, w5_6_fingers, w7_8_fingers

def size_ota_devices_from_kgm_and_currents(n_corner, p_corner, kgm_n, kgm_p, gain=50, bw=2e6):
    av= 50
    bw = 2e6
    gbw = bw * av
    therm_noise = 500e-9
    cload = 4e-12
    tan_thirty = math.tan(30*math.pi/180)
    alpha = 1/tan_thirty
    two_pi_alpha_gbw = 2*math.pi*alpha*gbw
    f2 = alpha*gbw
    total_current, m1_current, m6_current, beta_i_j, kcout_i_j, gain_i_j, thermal_rms_noise_i_j, beta_valid_i_j, gain_valid, thermal_noise_valid, kc_out = total_current_ota_v2(n_corner, p_corner, alpha, gbw, cload,
                                                                                                                           kgm_n, kgm_p, gain_spec=gain,
                                                                                                                           thermal_noise_spec=therm_noise)
    iden1_2 = n_corner.lookup(param1="kgm", param2="iden", param1_val=kgm_n)
    iden3_4 = p_corner.lookup(param1="kgm", param2="iden", param1_val=kgm_p)
    w1_2 = m1_current/iden1_2
    w3_4 = m1_current/iden3_4
    w5_6 = m6_current/iden1_2
    w7_8 = m6_current/iden3_4
    return w1_2, w3_4, w5_6, w7_8

if __name__ == "__main__":
    main()
