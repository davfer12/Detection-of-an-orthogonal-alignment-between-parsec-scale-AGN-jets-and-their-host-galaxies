from libraries import *


def linear(x,a,b):
    y = a*x+b
    return y    
    
def PA_difference(OPTICAL_PA,JET_PA):
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches = abs(OPTICAL_PA + 90 - JET_PA)
    index_90_180 = (90 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 180)
    index_180_270 = (180 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 270)
    index_270_360 = (270 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 360)
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_90_180] = 180 - X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_90_180]
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_180_270] = X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_180_270] - 180
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_270_360] = 360 - X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_270_360]
    return X_catalogue_Astrogeo_PA_DIFFERENCE_matches

def angle_conversion(angle):
    if not isinstance(angle,np.ndarray):
        angle = np.array([angle])
    angle_converted = np.fmod(angle+360*100,360)
    index_90_180 = (90 < angle_converted) & (angle_converted <= 180)
    index_180_270 = (180 < angle_converted) & (angle_converted <= 270)
    index_270_360 = (270 < angle_converted) & (angle_converted <= 360)
    angle_converted[index_90_180] = 180 - angle_converted[index_90_180]
    angle_converted[index_180_270] = angle_converted[index_180_270] - 180
    angle_converted[index_270_360] = 360 - angle_converted[index_270_360]
    return angle_converted  



def angle_conversion_major_PA(angle):
    if not isinstance(angle,np.ndarray):
        angle = np.array([angle])
    angle_converted = np.fmod(angle+360*100,180)
    angle_converted[angle_converted > 90] = angle_converted[angle_converted > 90] - 180
    return angle_converted
  
def mod(x_original,y_original,folded):
    dPA = y_original-x_original
    ii = dPA > 90
    jj = dPA < -90
    kk = dPA < 0
    x_mod = x_original.copy()
    y_mod = y_original.copy()
    y_mod[ii] -= 180
    y_mod[jj] += 180
    if folded:
        y_mod = x_mod + np.abs(y_mod-x_mod)
    return x_mod,y_mod
    
def MonteCarlo(angles,angles_err,N_draws,N_bins,N_bins_2):
    
    N_angles = len(angles)
    below_20 = np.where(angles_err < 20)[0]
    above_20 = np.where(angles_err > 20)[0]
    angles_rearranged = np.random.normal(loc=angles,scale=angles_err,size=(N_draws,N_angles))
    angles_rearranged = angle_conversion(angles_rearranged)

    bin_width = 90/N_bins


    bin_min = 0
    bin_max = 90
    Bins_array = np.arange(bin_min,bin_max+bin_width,bin_width)


    original_counts,original_bins=np.histogram(angles,bins=Bins_array)



    Histogram_counts = np.zeros((N_draws,N_bins))
    Histogram_below_20_counts = np.zeros((N_draws,N_bins))
    Histogram_above_20_counts = np.zeros((N_draws,N_bins))

    for i in range(N_draws):
        counts,bla=np.histogram(angles_rearranged[i],bins=Bins_array)
        counts_below_20,bla_below_20=np.histogram(angles_rearranged[i][below_20],bins=Bins_array,density=True)
        counts_above_20,bla_above_20=np.histogram(angles_rearranged[i][above_20],bins=Bins_array,density=True)
        Histogram_counts[i] = counts
        Histogram_below_20_counts[i] = counts_below_20
        Histogram_above_20_counts[i] = counts_above_20

    
    Cov_matrix_blue = np.cov(Histogram_counts,rowvar = False)
    
    Mean_counts = np.mean(Histogram_counts,axis=0)
    Mean_below_20_counts = np.mean(Histogram_below_20_counts,axis=0)
    Mean_above_20_counts = np.mean(Histogram_above_20_counts,axis=0)

    Errors_minus = np.zeros(N_bins)
    Errors_plus = np.zeros(N_bins)

    Errors_below_20_minus = np.zeros(N_bins)
    Errors_above_20_minus = np.zeros(N_bins)
    Errors_below_20_plus = np.zeros(N_bins)
    Errors_above_20_plus = np.zeros(N_bins)

    
    for i in range(N_bins):
        data = Histogram_counts[:,i]
        data_below_20 = Histogram_below_20_counts[:,i]
        data_above_20 = Histogram_above_20_counts[:,i]
        sigma1_left,median,sigma1_right = np.percentile(data,q=[15.8,50,100-15.8])
        sigma1_left_below_20,median_below_20,sigma1_right_below_20 = np.percentile(data_below_20,q=[15.8,50,100-15.8])
        sigma1_left_above_20,median_above_20,sigma1_right_above_20 = np.percentile(data_above_20,q=[15.8,50,100-15.8])
        Errors_minus[i] = Mean_counts[i] - sigma1_left
        Errors_below_20_minus[i] = Mean_below_20_counts[i] - sigma1_left_below_20
        Errors_above_20_minus[i] = Mean_above_20_counts[i] - sigma1_left_above_20
        Errors_plus[i] = sigma1_right - Mean_counts[i]
        Errors_below_20_plus[i] = Mean_below_20_counts[i] - sigma1_left_below_20
        Errors_above_20_plus[i] = Mean_above_20_counts[i] - sigma1_left_above_20
        average_sigma = np.mean([Errors_minus[i],Errors_plus[i]])
        average_sigma_below_20 = np.mean([Errors_below_20_minus[i],Errors_below_20_plus[i]])
        average_sigma_above_20 = np.mean([Errors_above_20_minus[i],Errors_above_20_plus[i]])
        x = np.linspace(np.min(data),np.max(data),100)
        curve = scipy.stats.norm.pdf(x,loc=Mean_counts[i],scale=average_sigma)

        x_below_20 = np.linspace(np.min(data_below_20),np.max(data_below_20),100)
        curve_below_20 = scipy.stats.norm.pdf(x_below_20,loc=Mean_below_20_counts[i],scale=average_sigma_below_20)

        x_above_20 = np.linspace(np.min(data_above_20),np.max(data_above_20),100)
        curve_above_20 = scipy.stats.norm.pdf(x_above_20,loc=Mean_above_20_counts[i],scale=average_sigma_above_20)    
    
    
    Bins_centers = (Bins_array[:-1] + Bins_array[1:])/2
    Y_errors = np.column_stack((Errors_minus,Errors_plus)).T

    return original_counts,Mean_counts,Y_errors,Bins_centers,Bins_array,Cov_matrix_blue    
    
def bin_N_significance(angles,angles_err,N_bins):
    original_counts,Mean_counts,Y_errors,Bins_centers,Bins_array,Cov_matrix_blue = MonteCarlo(angles,angles_err,10000,N_bins,20)
    errors = np.mean(Y_errors,axis=0)
    significance = (original_counts[1]-original_counts[0])/np.sqrt(np.sum(errors**2))
    return significance,errors

def cart2sph(x, y, z):
    xy = np.sqrt(x**2 + y**2)

    x_2 = x**2
    y_2 = y**2
    z_2 = z**2

    r = np.sqrt(x_2 + y_2 + z_2)
    theta = np.arctan2(xy,z)

    phi = np.arctan2(y,x) 

    return r, theta, phi

def R(u,v,theta):
    ux,uy,uz = u[0],u[1],u[2]
    vx,vy,vz = v[0],v[1],v[2]
    rot_matrix = np.array([[np.cos(theta)+ux**2*(1-np.cos(theta)),
                            ux*uy*(1-np.cos(theta))-uz*np.sin(theta),
                            ux*uz*(1-np.cos(theta))+uy*np.sin(theta)],
                           [uy*ux*(1-np.cos(theta))+uz*np.sin(theta),
                            np.cos(theta)+uy**2*(1-np.cos(theta)),
                            uy*uz*(1-np.cos(theta))-ux*np.sin(theta)],
                           [uz*ux*(1-np.cos(theta))-uy*np.sin(theta),
                            uz*uy*(1-np.cos(theta))+ux*np.sin(theta),
                           np.cos(theta)+uz**2*(1-np.cos(theta))]])
    
    v_new = np.matmul(rot_matrix,v)
    return v_new


def rot(spin,theta,azim):
    if len(spin.shape) == 1:
        spin = np.array([spin])
    inverse_cone = theta > 90
    azim = azim*np.pi/180
    theta = theta*np.pi/180
    spin_module = np.sqrt(np.sum(spin**2,axis=1))
     
    
    r = np.abs(np.tan(theta)*spin_module)
    a = spin[:,0]
    b = spin[:,1]
    c = spin[:,2]
    
    r_spin,theta_spin,phi_spin = cart2sph(a, b, c)
    
    
    new_theta = theta_spin+np.pi/2
    new_phi = phi_spin
    new_phi[new_theta > np.pi] = np.fmod(new_phi[new_theta > np.pi] + np.pi,2*np.pi)
    new_theta = np.fmod(theta_spin+np.pi/2,np.pi)

    
    perpx = r*np.cos(new_phi)*np.sin(new_theta)
    perpy = r*np.sin(new_phi)*np.sin(new_theta)
    perpz = -(a*perpx + b*perpy)/c
    perpz[c == 0] = 0
    
    
    
    rot_axis_x = spin[:,0]/spin_module
    rot_axis_y = spin[:,1]/spin_module
    rot_axis_z = spin[:,2]/spin_module
    
    rot_axis = np.array([rot_axis_x,rot_axis_y,rot_axis_z]).T
    perp = np.array([perpx,perpy,perpz]).T
    
    new_vector = np.zeros((len(spin),3))
    for i,(perp_i,rot_axis_i,azim_i) in enumerate(zip(perp,rot_axis,azim)):
        new_vector[i] = R(rot_axis_i,perp_i,azim_i)
    
    
    
    vector_final = spin + new_vector
    
    final_module = np.sqrt(np.sum(vector_final**2,axis=1))
    vector_final_x = (spin_module/final_module)*vector_final[:,0]
    vector_final_y = (spin_module/final_module)*vector_final[:,1]
    vector_final_z = (spin_module/final_module)*vector_final[:,2]
    
    vector_final = np.array([vector_final_x,vector_final_y,vector_final_z])
    vector_final = vector_final.T
    
    vector_final[inverse_cone] = -vector_final[inverse_cone]

    return vector_final