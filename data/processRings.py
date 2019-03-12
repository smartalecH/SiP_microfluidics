import numpy as np
import glob, os
from scipy.signal import find_peaks
from matplotlib import pyplot as plt


def characterize(wavelength,power,filename=""):
    ratio = 1/(wavelength[1] - wavelength[0])
    # find the peaks
    peaks, _ = find_peaks(power, height=0.01,width=30e-3*ratio,distance=2*ratio)

    z = np.polyfit(wavelength[peaks[1:]], np.diff(wavelength[peaks]), 2)
    p = np.poly1d(z)

    if filename is not "":
        plt.figure(figsize=(6,6))
        plt.subplot(2,1,1)
        plt.plot(wavelength,power)
        plt.plot(wavelength[peaks], power[peaks], "x")
        plt.grid(True)
        plt.ylabel('Power (au)')

        plt.subplot(2,1,2)
        plt.plot(wavelength,p(wavelength))
        plt.plot(wavelength[peaks[1:]],np.diff(wavelength[peaks]), "o")
        plt.grid(True)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('FSR (nm)')

        plt.tight_layout()
        plt.savefig(filename + ".png")
        plt.show()
    
    return p, peaks

# --------------------------------------------- #
# Get the baseline measurement
# --------------------------------------------- #
baselineFile = "baseline_drop.npz"
data = np.load(baselineFile)
wavelength = data['wavelength']
power = data['power']

pBaseline, pks = characterize(wavelength,power)




# --------------------------------------------- #
# Get all of the experimental measurements
# --------------------------------------------- #

#os.chdir("./")
experiments = []
concentration = [0]
FSR = [pBaseline(1560)]
peaks = [pks]
for file in glob.glob("*%*.npz"):
    experiments.append(file)
    concentration.append(file[0])
    data = np.load(file)
    wavelength = data['wavelength']
    power = data['power']
    prefix = file[0:-4]
    print(prefix)
    p, pks = characterize(wavelength,power,filename="")
    FSR.append(p(1560))
    peaks.append(pks)

peaks = np.array(peaks)

plt.figure(figsize=(6,5))
plt.plot(concentration,FSR,'o')
plt.grid(True)
plt.xlabel('Concentration NaCl (%)')
plt.ylabel('FSR (nm)')
plt.title('$\lambda_0=1560$ nm')
plt.tight_layout()
plt.savefig('fsrEvolution.png')




####################
wlFinal = []
concFinal = []
for k in range(peaks.size):
    temp = np.squeeze(np.array(peaks[k]))
    N = temp.size
    wlFinal = wlFinal + np.squeeze(wavelength[temp]).tolist()
    concFinal = concFinal + np.squeeze((np.tile(concentration[k],(N,)))).tolist()

concFinal = [float(line) for line in concFinal]
plt.figure(figsize=(6,5))
plt.plot(wlFinal,concFinal,'o')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Concentration NaCl (%)')
plt.tight_layout()
plt.grid(True)
plt.savefig('pksEvolution.png')
plt.show()
