import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Sample data (Gaussian distributed)
np.random.seed(0)
data = np.random.normal(loc=2, scale=1, size=1000)  # Example data with mean=2 and std=1

# Plot the histogram of the data
plt.hist(data, bins=30, density=True, alpha=0.5, color='b', label='Histogram')

# Compute the PDF values for a Gaussian distribution with the same mean and standard deviation
mean = np.mean(data)
std = np.std(data)
x = np.linspace(mean - 3*std, mean + 3*std, 1000)
pdf_values = norm.pdf(x, loc=mean, scale=std)

# Overlay the PDF curve on top of the histogram
plt.plot(x, pdf_values, color='r', label='Gaussian PDF')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Density')
plt.title('Histogram and Gaussian PDF')
plt.legend()

# Show plot
plt.show()