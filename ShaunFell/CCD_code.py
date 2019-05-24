import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import numpy as np
import glob
import astropy.constants as C
import astropy.units as U

bias_1 = "bias_1.FIT"
bias_2 = "bias_2.FIT"
bias_3 = "bias_3.FIT"


flats = glob.glob('Xavier/*_sec') + glob.glob('Xavier/*_sec.FIT')
containers = [fits.open(i) for i in flats]

mask = [700, 750, 1098, 1170]

flat_data = [np.asarray(i[0].data[mask[0]:mask[1], mask[2]:mask[3]], dtype='float64') for i in containers]

bias_3_container = fits.open(bias_3)
bias_2_container = fits.open(bias_2)
bias_1_container = fits.open(bias_1)

bias_1_data = np.asarray(bias_1_container[0].data[mask[0]:mask[1], mask[2]:mask[3]], dtype='float64')
bias_2_data = np.asarray(bias_2_container[0].data[mask[0]:mask[1], mask[2]:mask[3]], dtype='float64')
bias_3_data = np.asarray(bias_3_container[0].data[mask[0]:mask[1], mask[2]:mask[3]], dtype='float64')


#statistics of CCD

f1 = -9
f2 = -11

std_biases_difference = (bias_2_data - bias_1_data).std()
std_biases_difference_string = "The standard deviation of the difference of the two biases is {0}\n".format(std_biases_difference)
print(std_biases_difference_string)

std_flats_difference = (flat_data[f1] - flat_data[f2]).std()
std_flats_difference_string = "The standard deviation of the difference of two flats is {0}\n".format(std_flats_difference)
print(std_flats_difference_string)



mean_flats_addition = (flat_data[f1].mean() + flat_data[f2].mean())
mean_biases_addition = (bias_1_data.mean() + bias_2_data.mean())
mean_flats_addition_string = "The addition of the mean of two flats is {0}\n".format(mean_flats_addition)
mean_biases_addition_string = "The addition of the mean of two biases is {0}\n".format(mean_biases_addition)
print(mean_flats_addition_string)
print(mean_biases_addition_string)

Gain = (mean_flats_addition - mean_biases_addition)/(std_flats_difference**2 - std_biases_difference**2)
Gain_string = "The gain as calculated by equation 3 in the CCD Project packet is {0}".format(Gain)
print(Gain_string)

Read_noise = (Gain*std_biases_difference)/np.sqrt(2)
Read_noise_string = "The read noise as calculated by equation4 in the CCD Project packet is {0} in ADU".format(Read_noise)
Read_noise_string_e = "The read noise as calculated by equation4 in the CCD Project packet is {0} in electrons".format(Read_noise*Gain)
print(Read_noise_string)
print(Read_noise_string_e)




#linearity
means = np.asarray(list(map(np.mean, flat_data))) #Mean of each file
times = np.asarray([i[0].header['exptime'] for i in containers]) #Exposure time of each fits file
sort_mask = np.argsort(times) #Create sort mask
times = times[sort_mask] #apply sorting mask
means = means[sort_mask]
plt.title("Linearity")
plt.xlabel(r"Integration Time (s)")
plt.ylabel(r"counts (ADU's)")
plt.plot(times,means); plt.show() #Plot data


fig, ax = plt.subplots()
colLabels = ["Integration time (s)", "Mean Flat Level"]
cellText = [i for i in list(zip(times, means))]
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')
fig.tight_layout()
ax.table(cellText=cellText, colLabels = colLabels, loc='center'); plt.show()




#Dark Current


temp = np.asarray([-15, -10, -8, -6, -4.2, -2.3, 0, 2, 3.6, 5.8, 8.2, 10, 12.8, 16.2, 20]) + 273.15 #Kelvin +- 0.3 kelvin
Dark_Bias = np.asarray([13, 15, 16, 17, 19, 22, 24, 28, 32, 37, 43, 40, 62, 89, 139])*Gain#mean ADUs
gain = 2.3 #electrons per ADU
#k = C.k_B #Boltzmann constants in J/K
X = 1/temp
bias = 1000 #ADU
Dark_Current = bias + Dark_Bias
Dark_Current_func = lambda K, b, T: K*np.exp(-b*T)

plt.title("Dark Current versus Temperature")
plt.xlabel(r"$T^{-1}$ $(K^{-1})$")
plt.ylabel(r"Dark Current (electrons/pixel/sec)")
plt.plot(X, Dark_Current); plt.show()





