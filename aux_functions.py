from libraries import *


def linear(x, a, b):
    """
    Compute the value of a linear function.

    This function calculates the value of a linear function given by the formula:
        y = a * x + b

    Parameters:
    x (float or numpy.ndarray): The input value or array of values.
    a (float): The slope of the linear function.
    b (float): The y-intercept of the linear function.

    Returns:
    float or numpy.ndarray: The computed value(s) of the linear function.
    """
    y = a * x + b
    return y

    
def PA_difference(OPTICAL_PA, JET_PA):
    """
    Compute the absolute difference between optical and jet position angles.

    This function calculates the absolute difference between optical and jet position angles, taking into account the circular nature of angular measurements. The result is adjusted to ensure it falls within the range of 0 to 90 degrees.

    Parameters:
    OPTICAL_PA (numpy.ndarray or list): Array or list of optical position angles in degrees.
    JET_PA (numpy.ndarray or list): Array or list of jet position angles in degrees.

    Returns:
    numpy.ndarray: Array of absolute differences between optical and jet position angles, adjusted to the range [0, 180] degrees.

    Notes:
    - The difference is initially calculated as the absolute value of (OPTICAL_PA + 90 - JET_PA).
    - Adjustments are made for differences in the ranges [90, 180], [180, 270], and [270, 360] degrees to ensure the result is within [0, 90] degrees
    """
    # Calculate the initial absolute differences
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches = abs(OPTICAL_PA + 90 - JET_PA)
    
    # Adjust differences based on their range to ensure they are within [0, 180] degrees
    index_90_180 = (90 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 180)
    index_180_270 = (180 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 270)
    index_270_360 = (270 < X_catalogue_Astrogeo_PA_DIFFERENCE_matches) & (X_catalogue_Astrogeo_PA_DIFFERENCE_matches <= 360)
    
    # Apply adjustments to ensure differences are within [0, 180] degrees
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_90_180] = 180 - X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_90_180]
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_180_270] = X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_180_270] - 180
    X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_270_360] = 360 - X_catalogue_Astrogeo_PA_DIFFERENCE_matches[index_270_360]
    
    return X_catalogue_Astrogeo_PA_DIFFERENCE_matches

def process_coordinates(ra, dec):
        """
        Process right ascension and declination values: convert to radians and adjust
        for values greater than 180 degrees.
        """
        ra[ra > 180] -= 360
        dec[dec > 180] -= 360
        return ra * np.pi / 180, dec * np.pi / 180


def angle_conversion(angle):
    """
    Convert angle values to the range [-180, 180] degrees.

    This function adjusts angle values to ensure they fall within the range [-180, 180] degrees. The input angle can be a single value or a numpy array. The function first normalizes the angle to the range [0, 360) degrees, then adjusts it to ensure it is within [-180, 180] degrees.

    Parameters:
    angle (float, int, or numpy.ndarray): Angle value(s) in degrees. If a single value is provided, it will be converted to a numpy array.

    Returns:
    numpy.ndarray: Array of angles converted to the range [-180, 180] degrees.

    Notes:
    - The function normalizes angles to [0, 360) degrees first using modulo operation.
    - Angles in the ranges [90, 180], [180, 270], and [270, 360] are adjusted to fall within [-180, 180] degrees.
    """
    # Ensure input is a numpy array
    if not isinstance(angle, np.ndarray):
        angle = np.array([angle])
    
    # Normalize angles to the range [0, 360) degrees
    angle_converted = np.fmod(angle + 360 * 100, 360)
    
    # Identify and adjust angles based on their range to ensure they are within [0, 180] degrees
    index_90_180 = (90 < angle_converted) & (angle_converted <= 180)
    index_180_270 = (180 < angle_converted) & (angle_converted <= 270)
    index_270_360 = (270 < angle_converted) & (angle_converted <= 360)
    
    angle_converted[index_90_180] = 180 - angle_converted[index_90_180]
    angle_converted[index_180_270] = angle_converted[index_180_270] - 180
    angle_converted[index_270_360] = 360 - angle_converted[index_270_360]
    
    return angle_converted




def angle_conversion_major_PA(angle):
    """
    Convert angle values to the range [-90, 90] degrees.

    This function adjusts angle values to ensure they fall within the range [-90, 90] degrees. The input angle can be a single value or a numpy array. The function first normalizes the angle to the range [0, 180) degrees and then adjusts it to fall within [-90, 90] degrees.

    Parameters:
    angle (float, int, or numpy.ndarray): Angle value(s) in degrees. If a single value is provided, it will be converted to a numpy array.

    Returns:
    numpy.ndarray: Array of angles converted to the range [-90, 90] degrees.

    Notes:
    - The function normalizes angles to [0, 180) degrees first using modulo operation.
    - Angles greater than 90 degrees are adjusted to fall within the range [-90, 90] degrees.
    """
    # Ensure input is a numpy array
    if not isinstance(angle, np.ndarray):
        angle = np.array([angle])
    
    # Normalize angles to the range [0, 180) degrees
    angle_converted = np.fmod(angle + 360 * 100, 180)
    
    # Adjust angles greater than 90 degrees to fall within the range [-90, 90]
    angle_converted[angle_converted > 90] = angle_converted[angle_converted > 90] - 180
    
    return angle_converted
  
def mod(x_original, y_original, folded):
    """
    Adjust angles to ensure they are within a specified range and optionally fold the difference.

    This function adjusts the `y_original` values based on the difference between `x_original` and `y_original` angles. It handles cases where the difference exceeds 90 degrees by shifting the angles accordingly. If the `folded` parameter is `True`, it ensures that the resulting angles are always positive by folding the difference around `x_original`.

    Parameters:
    x_original (numpy.ndarray or list): Array or list of original x angles in degrees.
    y_original (numpy.ndarray or list): Array or list of original y angles in degrees.
    folded (bool): If `True`, the function will fold the difference to be positive. If `False`, it will return the adjusted angles without folding.

    Returns:
    numpy.ndarray: Adjusted x and y angles in degrees.
    """
    # Calculate the difference between y and x angles
    dPA = y_original - x_original
    
    # Determine where the difference exceeds ±90 degrees
    ii = dPA > 90
    jj = dPA < -90
    kk = dPA < 0
    
    # Create copies of original angles to modify
    x_mod = x_original.copy()
    y_mod = y_original.copy()
    
    # Adjust y_mod based on the calculated differences
    y_mod[ii] -= 180
    y_mod[jj] += 180
    
    # If folding is required, adjust y_mod to ensure it is always positive relative to x_mod
    if folded:
        y_mod = x_mod + np.abs(y_mod - x_mod)
    
    return x_mod, y_mod
    
def MonteCarlo(angles, angles_err, N_draws, N_bins):
    """
    Perform a Monte Carlo simulation to analyze the distribution of angles and their errors.

    This function generates multiple realizations of angle distributions based on given angles and their associated errors. It calculates histograms for the simulated distributions and computes statistical measures such as means and errors. The function also generates Gaussian curves based on these statistics.

    Parameters:
    angles (numpy.ndarray): Array of angle values in degrees.
    angles_err (numpy.ndarray): Array of angle error values in degrees.
    N_draws (int): Number of Monte Carlo draws.
    N_bins (int): Number of bins for the histogram.

    Returns:
    tuple: Contains the following elements:
        - original_counts (numpy.ndarray): Histogram counts of the original angles.
        - Mean_counts (numpy.ndarray): Mean histogram counts from the Monte Carlo simulations.
        - Y_errors (numpy.ndarray): Array of errors for each bin.
        - Bins_centers (numpy.ndarray): Centers of the histogram bins.
        - Bins_array (numpy.ndarray): Array defining the histogram bins.
        - Cov_matrix_blue (numpy.ndarray): Covariance matrix of the histogram counts.
    """

    # Determine the number of angles
    N_angles = len(angles)
        
    # Generate Monte Carlo simulations
    angles_rearranged = np.random.normal(loc=angles, scale=angles_err, size=(N_draws, N_angles))
    angles_rearranged = angle_conversion(angles_rearranged)

    # Define histogram bin edges
    bin_width = 90 / N_bins
    bin_min = 0
    bin_max = 90
    Bins_array = np.arange(bin_min, bin_max + bin_width, bin_width)

    # Calculate histogram for the original angles
    original_counts, original_bins = np.histogram(angles, bins=Bins_array)

    # Initialize arrays for histogram counts from simulations
    Histogram_counts = np.zeros((N_draws, N_bins))

    # Compute histograms for each draw
    for i in range(N_draws):
        counts, _ = np.histogram(angles_rearranged[i], bins=Bins_array)
        Histogram_counts[i] = counts


    # Compute the covariance matrix of the histogram counts
    Cov_matrix_blue = np.cov(Histogram_counts, rowvar=False)
    
    # Calculate mean histogram counts
    Mean_counts = np.mean(Histogram_counts, axis=0)


    # Initialize error arrays
    Errors_minus = np.zeros(N_bins)
    Errors_plus = np.zeros(N_bins)


    # Calculate errors for each bin
    for i in range(N_bins):
        data = Histogram_counts[:, i]
     
        sigma1_left, median, sigma1_right = np.percentile(data, q=[15.8, 50, 100 - 15.8])

        
        Errors_minus[i] = Mean_counts[i] - sigma1_left
        Errors_plus[i] = sigma1_right - Mean_counts[i]

        
        # Average error calculations
        average_sigma = np.mean([Errors_minus[i], Errors_plus[i]])
  

    # Calculate bin centers
    Bins_centers = (Bins_array[:-1] + Bins_array[1:]) / 2

    # Stack errors for return
    Y_errors = np.column_stack((Errors_minus, Errors_plus)).T

    return original_counts, Mean_counts, Y_errors, Bins_centers, Bins_array, Cov_matrix_blue
    
    
def bin_N_significance(angles, angles_err, N_bins):
    """
    Calculate the significance of the difference between histogram bins based on Monte Carlo simulations.

    This function performs a Monte Carlo simulation to generate histograms of angles and calculates the significance of the difference between adjacent bins.

    Parameters:
    angles (numpy.ndarray): Array of angle values in degrees.
    angles_err (numpy.ndarray): Array of angle error values in degrees.
    N_bins (int): Number of bins for the histogram.

    Returns:
    tuple: Contains the following elements:
        - significance (float): Significance of the difference between the first two bins.
        - errors (numpy.ndarray): Array of errors for each bin.
    """

    # Perform Monte Carlo simulations
    original_counts, Mean_counts, Y_errors, Bins_centers, Bins_array, Cov_matrix_blue = MonteCarlo(angles, angles_err, 10000, N_bins)

    # Calculate errors as the mean of the plus and minus errors
    errors = np.mean(Y_errors, axis=0)

    # Calculate significance of the difference between the first two bins
    significance = (original_counts[1] - original_counts[0]) / np.sqrt(np.sum(errors**2))

    return significance, errors

def cart2sph(x, y, z):
    """
    Convert Cartesian coordinates to spherical coordinates.

    Parameters:
    x (numpy.ndarray or float): Cartesian x-coordinate(s).
    y (numpy.ndarray or float): Cartesian y-coordinate(s).
    z (numpy.ndarray or float): Cartesian z-coordinate(s).

    Returns:
    tuple: (r, theta, phi)
        - r (numpy.ndarray or float): Radial distance from the origin.
        - theta (numpy.ndarray or float): Polar angle in radians (angle from the z-axis).
        - phi (numpy.ndarray or float): Azimuthal angle in radians (angle in the xy-plane from the x-axis).
    """
    
    xy = np.sqrt(x**2 + y**2)

    # Radial distance
    r = np.sqrt(x**2 + y**2 + z**2)
    
    # Polar angle (angle from the z-axis)
    theta = np.arctan2(xy, z)
    
    # Azimuthal angle (angle in the xy-plane from the x-axis)
    phi = np.arctan2(y, x)
    
    return r, theta, phi

def R(u, v, theta):
    """
    Rotate vector v around axis u by angle theta.

    Parameters:
    u (numpy.ndarray): 3D vector representing the axis of rotation (should be a unit vector).
    v (numpy.ndarray): 3D vector to be rotated.
    theta (float): Angle of rotation in radians.

    Returns:
    numpy.ndarray: The rotated vector.
    """
    
    # Extract components of u and v
    ux, uy, uz = u[0], u[1], u[2]
    vx, vy, vz = v[0], v[1], v[2]
    
    # Construct the rotation matrix using Rodrigues' rotation formula
    rot_matrix = np.array([
        [np.cos(theta) + ux**2 * (1 - np.cos(theta)),
         ux*uy * (1 - np.cos(theta)) - uz * np.sin(theta),
         ux*uz * (1 - np.cos(theta)) + uy * np.sin(theta)],
        [uy*ux * (1 - np.cos(theta)) + uz * np.sin(theta),
         np.cos(theta) + uy**2 * (1 - np.cos(theta)),
         uy*uz * (1 - np.cos(theta)) - ux * np.sin(theta)],
        [uz*ux * (1 - np.cos(theta)) - uy * np.sin(theta),
         uz*uy * (1 - np.cos(theta)) + ux * np.sin(theta),
         np.cos(theta) + uz**2 * (1 - np.cos(theta))]
    ])
    
    # Apply the rotation matrix to vector v
    v_new = np.matmul(rot_matrix, v)
    
    return v_new


def rot(spin, theta, azim):
    """
    Rotate 3D vectors according to specified angles.

    Parameters:
    spin (numpy.ndarray): 2D array where each row is a 3D vector to be rotated.
    theta (float or numpy.ndarray): Angle of rotation in degrees.
    azim (float or numpy.ndarray): Azimuthal angle in degrees.

    Returns:
    numpy.ndarray: Rotated 3D vectors.
    """
    if len(spin.shape) == 1:
        spin = np.array([spin])
        
    # Convert angles from degrees to radians
    theta = np.radians(theta)
    azim = np.radians(azim)

    # Identify vectors with theta > 90
    inverse_cone = theta > np.pi / 2

    # Vector magnitudes
    spin_module = np.linalg.norm(spin, axis=1)
    
    # Compute spherical coordinates for spin vectors
    a, b, c = spin[:,0], spin[:,1], spin[:,2]
    r, theta_spin, phi_spin = cart2sph(a, b, c)
    
    # Adjust theta and phi for rotation
    new_theta = np.fmod(theta_spin + np.pi / 2, np.pi)
    new_phi = np.fmod(phi_spin + np.where(new_theta > np.pi, np.pi, 0), 2 * np.pi)
    
    # Compute perpendicular vector components
    r_spin = np.abs(np.tan(theta) * spin_module)
    perpx = r_spin * np.cos(new_phi) * np.sin(new_theta)
    perpy = r_spin * np.sin(new_phi) * np.sin(new_theta)
    perpz = np.where(c != 0, -(a * perpx + b * perpy) / c, 0)
    
    # Rotation axis (normalized)
    rot_axis = spin / spin_module[:, np.newaxis]
    
    # Perform rotation
    new_vector = np.array([R(axis, [px, py, pz], az) for axis, px, py, pz, az in zip(rot_axis, perpx, perpy, perpz, azim)])
    
    # Adjust the final vector's magnitude and direction
    final_vector = spin + new_vector
    final_module = np.linalg.norm(final_vector, axis=1)
    scale_factor = spin_module / final_module
    vector_final = final_vector * scale_factor[:, np.newaxis]
    
    # Handle vectors initially outside the cone
    vector_final[inverse_cone] = -vector_final[inverse_cone]

    return vector_final